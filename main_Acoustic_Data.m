%% ========================================================================
%  Acoustic Synthetic Data Generation Using Modified Nearly Analytical
%  Discritization (MNAD) Numerical Solver
% ========================================================================
%  Purpose:
%    Generates synthetic acoustic seismic data. Loads a velocity
%    model, builds a PML-padded domain, sets up source/receiver
%    acquisition geometry, builds the MNAD coefficient
%    matrices, and forward-models shot gathers in parallel.
%
%  Compatible with both synthetic acquisitions OBS / streamer geometries
%
%  Outputs:
%    - Illuminations             : output directories icluding illumination 
%    - Wavefields                : output directories icluding nt/nf/2 time
%      slice of wavefield
%    - data_Acoustic.mat         : forward-modeled shot gather data
%    - profile                   : key parameters into a profile struct 
% ========================================================================
%  Author: Toktam Zand
%  Contact: tzand@igf.edu.pl

clear ; clc ; %close all ;

IW = true ; % Create output folders for illumination maps and wavefield
            % snapshots (used by later steps in the RTM workflow, 
            % not by this script directly). If you no need change IW =
            % false ;
if IW == true
    if ~exist('Illuminations','dir'), mkdir Illuminations; end
    if ~exist('Wavefields','dir'), mkdir Wavefields; end
end

% % Set the number inside parpool based on your system
poo = gcp('nocreate');
if isempty(poo), parpool(2);clear poo; end

% Log start time for run-time tracking
date = string(datetime('now')) ;
formatSpec = 'Start -> Date: %s \n';
fprintf(formatSpec,date)

% % Thickness of PML (Perfectly Matched Layer) absorbing boundary, in grid cells
d_pml = 40 ;

% % ______________________________________________________________________% %
% %  ---------------------- Synthetic Data Generation ------------------

%% --- Load velocity model ---
newname = struct2cell(load('Velocity.mat')) ;
V = newname{1} ;
[nz,nx] = size(V) ;            % Number of grid points in x and z directions
V_min = min(V(:)) ;
V_max = max(V(:)) ;

%% --- Frequency content ---
fc = 10 ;                       % Dominant frequency of the Ricker wavelet [Hz]
f_max  = 3.*fc ;               % Approximate maximum frequency content (for e.g. stability/dispersion checks

profile.fc = fc ;    
profile.f_max = f_max ; 

%% --- Domain Discritization ---
% Compute stable time step [s] and spatial steps [m]

dx = V_min/(2.5*f_max) ;
dz = dx ;
% % dx and dz may change regarding your need but must be smaller than what
% calculated above corresponding maximum velocity and frequency 

dt = (0.65*min(dx,dz))/V_max ;


%% --- Time axis ---
% The recording time
T_max = 2 ;                               % Set the total recording duration in seconds
dt = 0.8.*dt ;                            % Time step [s] 
t = 0:dt:T_max-dt ;                       % Time vector
nt = length(t) ;                          % Number of time samples
fs = 1/dt ;                               % Sampling rate [Hz]

profile.dt = dt ;        
profile.T_max = T_max ;        

%% --- Spatial grid (x-direction) ---
dx = floor(dx) ;                         % Grid spacing in x [m]
x = (0:nx-1)*dx ;

profile.dx = dx ;         
profile.nx = nx ;         

%% --- Spatial grid (z-direction) ---
dz = floor(dz) ;                          % Number of grid points and spacing [m] in z
z = (0:nz-1)*dz ;

profile.dz = dz ;         
profile.nz = nz ;        

%% --- Acquisition geometry (synthetic) ---
zsc = 2*dz ;                                   % Set sources depth [m]
zrc = 5*dz ;                                   % Set receivers depth [m]
h_offset = 2000 ;                              % Set the horizontal offset of receiver array [m] if needed
ds = 20*dx ;                                    % Source spacing in lateral direction [m]
dr = 1*dx ;                                   % Receiver spacing in lateral direction [m]
slag = 2*dx ;                               % Source lag from the sides of model in lateral direction [m]
rlag = dx ;                                 % Source lag from the sides of model in lateral direction [m]

% Build source/receiver coordinates (in meters) for all shots
[Source, Receiver] = acquisi(V,dx,dz,zsc,zrc,h_offset,ds,dr,slag,rlag) ;
ns = length(Source) ;                                   % The number of source shots

% Preallocate per-shot coordinate cells
xs = cell(1,ns) ;
zs = cell(1,ns) ;
xr = cell(1,ns) ;
zr = cell(1,ns) ;

% Convert source/receiver positions from meters to grid indices,
% and shift by d_pml to account for the PML padding added to the model
for ss = 1:ns
    xs{ss} = d_pml + Source(1,ss)./dx ;        % Source x grid index
    zs{ss} = d_pml + Source(2,ss)./dz ;        % Source z grid index

    xr{ss} = d_pml + Receiver{1,ss}./dx ;      % Receiver x grid indices
    zr{ss} = d_pml + Receiver{2,ss}./dz ;      % Receiver z grid indices
end

profile.ns = ns ;
profile.Source = Source ;
profile.Receiver = Receiver ;

%% --- Source wavelet ---
% Generate Ricker wavelet
t0 = 1.1/fc ;        % [s]   delay wavelet
tau = pi*(t-t0)*fc ;         % [-]   Ricker wavelet argument
w  = ((1 - 2*tau.^2) .* exp( -tau.^2 ))' ; % Analytical Ricker Wavelet solution

profile.w = w ;

%% --- Pad velocity model with PML boundary layers ---
% Extend the model edges outward by d_pml cells on each side (x then z),
% by repeating the edge values — this is a common way to build the PML
% "buffer" region  around the physical domain.
v = [repmat(V(:,1),1,d_pml) , V , repmat(V(:,end),1,d_pml)] ;
v = [repmat(v(1,:),d_pml,1) ; v ; repmat(v(end,:),d_pml,1)] ;

%% --- Build MNAD coefficient matrices ---
% Each *_gen function builds one element of the MNAD kernal tensor A used to
% propagate acoustic wavefield "u" and gradient fields vs x and z axces "v" and "w", respectively. 
% The cofficients provided in *_gen can be find in Appendix A of the
% manuscript

Auu = Auu_gen(dx,dz,dt,v) ;
Auv = Auv_gen(dx,dz,dt,v) ;
Auw = Auw_gen(dx,dz,dt,v) ;
Avu = Avu_gen(dx,dz,dt,v) ;
Avv = Avv_gen(dx,dz,dt,v) ;
Avw = Avw_gen(dx,dz,dt,v) ;
Awu = Awu_gen(dx,dz,dt,v) ;
Awv = Awv_gen(dx,dz,dt,v) ;
Aww = Aww_gen(dx,dz,dt,v) ;
save As Auu Auv Auw Avu Avv Avw Awu Awv Aww

%% --- Forward model shot gathers in parallel ---
% Each worker forward-propagates one shot and returns its recorded data
% Function "fwd_PDE_exact_MNAD" solve the Equation 9 in the manuscript.

data = cell(1,ns);
parfor ss = 1:ns
    data{ss} = fwd_PDE_exact_MNAD(dx,xs{ss},xr{ss},dz,zs{ss},zr{ss},v,dt,nt,w,d_pml,ss,Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww) ;
end

% Log completion time
date = string(datetime('now')) ;
formatSpec = 'Data generated -> Date: %s \n';
fprintf(formatSpec,date)


save data_MNAD_Acoustic data profile
