function data = fwd_PDE_exact_MNAD(dx,xs,xr,dz,zs,zr,v,dt,nt,w,d_pml,ss,Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww,indx)
% % ========================================================================
%  Acoustic Wave Propagation Solver (MNAD scheme)
% ========================================================================
%  Purpose:
%    Propagates a single shot's wavefield forward in time using the
%    Modified Nearly Analytical Discretization (MNAD) scheme, injects a
%    point source, records receiver data, and (optionally) saves
%    snapshots of the wavefield's second time-derivative and its
%    illumination map for later use in RTM/LSRTM gradient computation.
%
%  Inputs:
%    dx, dz     - Grid spacing in x and z [m]
%    xs, zs     - Source grid index (x, z)
%    xr, zr     - Receiver grid indices (x, z)
%    v          - Padded (PML-included) velocity model
%    dt         - Time step [s]
%    nt         - Number of time samples
%    w          - Source wavelet, length nt
%    d_pml      - PML thickness in grid cells 
%    ss         - Shot index 
%    indx       - If 1, save per-shot second time-derivative wavefield snapshots (Wavefields/)
%                 and illumination map (Illuminations/). If 0/omitted,
%                 skip both (saves time/disk for shots where this isn't
%                 needed, e.g. shots not used for gradient computation).
%                 NOTE: defaults to 0 if not supplied 
%
%  Outputs:
%    data       - {1x3} cell: {datau, datav, dataw}, each [nt x nr],
%                 recorded receiver data for the pressure (u), vertical gradient of the preesure (v), 
%                 and horizontal gradient of the preesure (w) wavefields
%                 components respectively
% ========================================================================

if nargin < 22 || isempty(indx)
    indx = 0 ;   % Default: do not save wavefields/illumination
end

V2dt2 = v.^2.*(dt^2) ;


% % --------- Snapshot storage sizing 
% Number of wavefield snapshots to retain based on required
% sampling rate by imaging condition, scaled down as the
% max/min velocity ratio grows (need to control memory use)
nf = floor(nt/(2*3*max(v(:))/min(v(:)))) ;
T = floor(nt/nf/2);

% Preallocate storage for the second-time-derivative wavefield snapshots
% (interior region only, i.e. PML stripped out) — used later for the
% RTM/LSRTM imaging condition
[nz,nx] = size(v) ;
UU = zeros((nz-2*d_pml)*(nx-2*d_pml),2*nf + 1) ;
ttt = 1 ;
illum = 0 ;    % Source-side illumination approximation (sum of squared second time-derivative of wavefield)

% % --------- Receiver data storage 
datau = zeros(nt,length(xr)) ;
datav = zeros(nt,length(xr)) ;
dataw = zeros(nt,length(xr)) ;

% % --------- Wavefield state variables 
% U, V, W: current wavefield components; pU, pV, pW: previous time step
pU = zeros(nz,nx) ; U = zeros(nz,nx) ;
pV = zeros(nz,nx) ; V = zeros(nz,nx) ;
pW = zeros(nz,nx) ; W = zeros(nz,nx) ;

% % --------- PML auxiliary variables 
% MNAD/CPML absorbing boundaries require auxiliary "memory" fields per
% wavefield component (U, V, W), per direction (x, z), and per boundary
% side (left/right for x, top/bottom for z).
U1_xl = zeros(nz,nx) ; U2_xl = zeros(nz,nx) ; U3_xl = zeros(nz,nx) ;
U1_xr = zeros(nz,nx) ; U2_xr = zeros(nz,nx) ; U3_xr = zeros(nz,nx) ;
U1_zt = zeros(nz,nx) ; U2_zt = zeros(nz,nx) ; U3_zt = zeros(nz,nx) ;
U1_zb = zeros(nz,nx) ; U2_zb = zeros(nz,nx) ; U3_zb = zeros(nz,nx) ;

V1_xl = zeros(nz,nx) ; V2_xl = zeros(nz,nx) ; V3_xl = zeros(nz,nx) ;
V1_xr = zeros(nz,nx) ; V2_xr = zeros(nz,nx) ; V3_xr = zeros(nz,nx) ;
V1_zt = zeros(nz,nx) ; V2_zt = zeros(nz,nx) ; V3_zt = zeros(nz,nx) ;
V1_zb = zeros(nz,nx) ; V2_zb = zeros(nz,nx) ; V3_zb = zeros(nz,nx) ;

W1_xl = zeros(nz,nx) ; W2_xl = zeros(nz,nx) ; W3_xl = zeros(nz,nx) ;
W1_xr = zeros(nz,nx) ; W2_xr = zeros(nz,nx) ; W3_xr = zeros(nz,nx) ;
W1_zt = zeros(nz,nx) ; W2_zt = zeros(nz,nx) ; W3_zt = zeros(nz,nx) ;
W1_zb = zeros(nz,nx) ; W2_zb = zeros(nz,nx) ; W3_zb = zeros(nz,nx) ;

% % --------- PML damping coefficients 
R = 1*10^-4 ;     % Target PML reflection coefficient
[c1xl,c2xl,c3xl,c1xr,c2xr,c3xr,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb] = C_cof(R,d_pml,dx,dz,dt,nx,nz,v) ;

% % --------- Time-stepping loop 
for tt = 1:nt
    
    %%%%%%%%%%%%%%% U: compute PML-corrected spatial derivatives
    [U1_x,U1_xl,U2_xl,U3_xl,U1_xr,U2_xr,U3_xr] =...
        Ux123_gen(U1_xl,U2_xl,U3_xl,U1_xr,U2_xr,U3_xr,U,dx,nz,nx,c1xl,c2xl,c3xl,c1xr,c2xr,c3xr) ;
    [U1_z,U1_zt,U2_zt,U3_zt,U1_zb,U2_zb,U3_zb] =...
        Uz123_gen(U1_zt,U2_zt,U3_zt,U1_zb,U2_zb,U3_zb,U,dz,nz,nx,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb) ;
    %%%%%%%%%%%%%%% V: compute PML-corrected spatial derivatives
    [V1_x,V1_xl,V2_xl,V3_xl,V1_xr,V2_xr,V3_xr] =...
        Ux123_gen(V1_xl,V2_xl,V3_xl,V1_xr,V2_xr,V3_xr,V,dx,nz,nx,c1xl,c2xl,c3xl,c1xr,c2xr,c3xr) ;
    [V1_z,V1_zt,V2_zt,V3_zt,V1_zb,V2_zb,V3_zb] =...
        Uz123_gen(V1_zt,V2_zt,V3_zt,V1_zb,V2_zb,V3_zb,V,dz,nz,nx,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb) ;
    %%%%%%%%%%%%%%% W: compute PML-corrected spatial derivatives
    [W1_x,W1_xl,W2_xl,W3_xl,W1_xr,W2_xr,W3_xr] =...
        Ux123_gen(W1_xl,W2_xl,W3_xl,W1_xr,W2_xr,W3_xr,W,dx,nz,nx,c1xl,c2xl,c3xl,c1xr,c2xr,c3xr) ;
    [W1_z,W1_zt,W2_zt,W3_zt,W1_zb,W2_zb,W3_zb] =...
        Uz123_gen(W1_zt,W2_zt,W3_zt,W1_zb,W2_zb,W3_zb,W,dz,nz,nx,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb) ;
    
    % --------- Inject point source into U at (zs, xs) at this time step
    U(zs,xs) = U(zs,xs) + V2dt2(zs,xs)*w(tt) ;
    
    % --------- MNAD update equations: advance U, V, W one time step
    nU = reshape(Auu*U(:),nz,nx) + reshape(Auv*V(:),nz,nx) + reshape(Auw*W(:),nz,nx) - pU - V2dt2.*(U1_x + U1_z) ;
    nV = reshape(Avv*V(:),nz,nx) + reshape(Avu*U(:),nz,nx) + reshape(Avw*W(:),nz,nx) - pV - V2dt2.*(V1_x + V1_z) ;
    nW = reshape(Aww*W(:),nz,nx) + reshape(Awu*U(:),nz,nx) + reshape(Awv*V(:),nz,nx) - pW - V2dt2.*(W1_x + W1_z) ;

    % --------- Record receiver data at this time step
    for rr = 1:length(xr)
        datau(tt,rr) = U(zr(rr),xr(rr)) ;
        datav(tt,rr) = V(zr(rr),xr(rr)) ;
        dataw(tt,rr) = W(zr(rr),xr(rr)) ;
    end
    
    % --------- Optionally save wavefield snapshots + illumination
    if indx == 1
        
        if tt==1
            d2u = nU/dt^2;
        elseif tt==2
            d2u = (nU - 2*U)/dt^2;
        else
            d2u = (nU - 2*U + pU)/dt^2;
        end
        
        % --------- Keep the wavefield every T steps
        if rem(tt,T) == 0
            internal = d2u(d_pml+1:end-d_pml,d_pml+1:end-d_pml) ;
            UU(:,ttt) = internal(:) ;
            illum = illum + internal.^2 ;
            ttt = ttt + 1 ;
        end
    end
    
    % --------- Advance time level:
    pU = U ;
    U = nU ;
    
    pV = V ;
    V = nV ;
    
    pW = W ;
    W = nW ;
    
         if mod(tt,50)==0, imagesc(U);caxis(0.1.*caxis);drawnow;end

end

% % --------- Save per-shot outputs (only if indx == 1)
if indx == 1
    save(['./Wavefields/UU_' num2str(ss)],'UU')
    illum = T*illum ;
    save(['./Illuminations/illum_' num2str(ss)],'illum')
end

data{1} = datau ;
% data{2} = datav ;   % Use this if you need the lateral gradient of pressure field to be recorded 
% data{3} = dataw ;   % Use this if you need the vertical gradient of pressure field to be recorded

end


