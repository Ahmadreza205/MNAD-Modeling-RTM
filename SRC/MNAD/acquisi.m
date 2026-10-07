function [Source, Receiver] = acquisi(V,dx,dz,zsc,zrc,h_offset,ds,dr,slag,rlag,acqType,doSave,manualGeom)
%% ========================================================================
%  Acquisition Geometry Generator
% ========================================================================
%  Purpose:
%    Builds source and receiver coordinates (in meters) for a synthetic
%    seismic acquisition, given a velocity model and geometry parameters.
%    Sources are placed at evenly spaced x-locations across the model;
%    receiver layout depends on the chosen acquisition type.
%
%  Inputs:
%    V           - Velocity model, size [nz x nx]
%    dx, dz      - Grid spacing in x and z [m]
%    zsc, zrc    - Source and receiver depth [m] (converted to grid index
%                  internally; unused in 'arbitrary' mode)
%    h_offset    - Horizontal receiver offset [m] (used by streamer modes)
%    ds, dr      - Source and Receiver spacing in lateral direction [m]
%    slag, rlag  - Source and Receiver lags from the sides of model in lateral direction [m]
%    acqType     - (optional) Acquisition geometry type, one of:
%                  'fixspread' (default)   - Full-fold split-spread streamer,
%                                             every source sees all nx receivers
%                  'onesided'              - One-sided streamer, receivers
%                                             span [xs, xs+h_offset]
%                  'splitsided'            - Two-sided streamer, receivers
%                                             span [xs-h_offset, xs+h_offset]
%                  'obn'                   - Ocean-bottom-node geometry
%                  'arbitrary'             - Manually specified geometry,
%                                             supplied via manualGeom
%    doSave    - (optional) If true, saves Source/Receiver to .mat files
%                in the current directory. Default: false.
%    manualGeom- (required only if acqType='arbitrary') Struct with grid
%                index fields:
%                  .xs - [1 x ns] source x grid indices
%                  .zs - [1 x ns] source z grid indices
%                  .xr - {1 x ns} cell, receiver x grid indices per shot
%                  .zr - {1 x ns} cell, receiver z grid indices per shot
%                All indices are in the SAME grid-index convention as the
%                other modes (i.e. not yet converted to meters — this
%                function still does that conversion at the end).
%
%  Outputs:
%    Source    - [2 x ns] matrix of [x;z] source coordinates in meters
%    Receiver  - {2 x ns} cell array of [x;z] receiver coordinates in
%                meters, one column pair per shot
% ========================================================================

if nargin < 11 || isempty(acqType), acqType = 'fixspread' ;    end   
if nargin < 12 || isempty(doSave),  doSave  = false ;          end
if nargin < 13,                     manualGeom = [] ;          end

date = string(datetime('now')) ;
fprintf('Acquisition-> Date: %s \n', date)

[nz,nx] = size(V) ; 
ds = fix(ds/dx) ;
dr = fix(dr/dx) ;

if ds == 0, ds = 1 ; end 
if dr == 0, dr = 1 ; end

switch lower(acqType)

    case 'arbitrary'
        %% Manually specified geometry — bypasses automatic layout entirely
        if isempty(manualGeom)
            error('acquisi:missingManualGeom', ...
                  'acqType is ''arbitrary'' but manualGeom was not provided.') ;
        end
        req = {'xs','zs','xr','zr'} ;
        for kk = 1:numel(req)
            if ~isfield(manualGeom,req{kk})
                error('acquisi:missingManualField', ...
                      'manualGeom is missing required field ''%s''.', req{kk}) ;
            end
        end

        xs = manualGeom.xs ;
        zs = manualGeom.zs ;
        xr = manualGeom.xr ;
        zr = manualGeom.zr ;

    otherwise
        % Source x-locations: evenly spaced across the model width.
        % To use custom (manually chosen) source positions instead,
        % either edit the line below directly, or use acqType =
        % 'arbitrary' with manualGeom for a cleaner, non-hardcoded
        % way to pass in custom positions.   
        
        xs = ceil(slag./dx):ds:(nx-floor(slag./dx)) ;   
        zs = zeros(size(xs)) + floor(zsc./dz) ;
        
        ns = length(xs) ;
        
        xr = cell(1,ns) ;
        zr = cell(1,ns) ;
        h_off = floor(h_offset./dx) ;

        switch lower(acqType)

            case 'onesided'
                %% One-sided streamer: receivers extend from the source out to h_offset
                for ss = 1:ns
                    if xs(ss) + h_off  > nx - floor(rlag./dx)
                        xr{ss} = xs(ss):dr:(nx-floor(rlag./dx)) ;
                    else
                        xr{ss} = xs(ss):dr:xs(ss)+h_off-1 ;
                    end
                    zr{ss} = 0.*xr{ss} + floor(zrc./dz) ;
                end
                
            case 'splitspread'
                %% Split-spread streamer: receivers extend from the source out to h_offset from both sides
                for ss = 1:ns
                    if xs(ss) + h_off  > nx - floor(rlag./dx)
                        xr{ss} = (xs(ss) - h_off):dr:(nx - floor(rlag./dx)) ;
                    elseif xs(ss) - h_off  <= 0
                        xr{ss} = 1:dr:(xs(ss)+h_off-1) ;
                    else
                        xr{ss} = (xs(ss) - h_off):dr:(xs(ss)+h_off-1) ;
                    end
                    zr{ss} = 0.*xr{ss} + floor(zrc./dz) ;
                end

            case 'obn'
                %% Ocean-Bottom-Node geometry 
                %  NOTE: This mode originally depended on V_water
                %  (water velocity, used to detect the seafloor depth per column).
                %  V_water consider as the value of the first grid of
                %  velocity V

                V_water = V(1,1) ;
                clear xz zs
                xs = ceil(rlag./dx):dr:(nx - floor(rlag./dx)) ;
                ns = length(xs) ;
                for rr = 1:ns
                    col_V = V(:,xs(rr)) ;
                    indx_OBN = find(col_V==V_water) ;
                    zs(rr) = max(indx_OBN) ;
                end
                
                for ss = 1:ns
                    xr{ss} = ceil(slag./dx):ds:(nx - floor(slag./dx)) ;
                    zr{ss} = 0.*xr{ss} + floor(zsc./dz) ;
                end

            case 'fixspread'
                %% Fix-spread streamer: every source records at every receiver
                %  position across the full model width (full-fold geometry)
                for ss = 1:ns
                    xr{ss} = ceil(rlag./dx):dr:(nx-floor(rlag./dx)) ;
                    zr{ss} = 0.*xr{ss} + floor(zrc./dz) ;
                end
            otherwise
                error('acquisi:unknownType', ...
                      ['Unknown acqType "%s". Use ''splitspread'', ', ...
                       '''onesided'', ''obn'', or ''arbitrary''.'], acqType) ;
        end
end

%% --- Convert grid indices to meters ---
Source = [dx.*xs ; dz.*zs] ;

xr_l = cell(1,ns) ;
zr_l = cell(1,ns) ;
for ss = 1:ns
    xr_l{ss} = dx.*xr{ss} ;
    zr_l{ss} = dz.*zr{ss} ;
end
Receiver = [xr_l ; zr_l] ;

if doSave
    save Source Source
    save Receiver Receiver
end

end
