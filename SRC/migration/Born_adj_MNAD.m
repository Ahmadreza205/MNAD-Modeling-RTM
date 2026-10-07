function [Im]= Born_adj_MNAD(dx,dz,nt,dt,v,xr,zr,data,d_pml,ss,Auu,Auv,Auw,Avu,Avv,Avw,Awu,Awv,Aww)
% % ========================================================================
%  Born Adjoint Operator (RTM Migration Kernel, MNAD scheme)
% ========================================================================
%  Purpose:
%    Back-propagates the scattered-data residual from receiver locations
%    backward in time through the background velocity model, then
%    cross-correlates the resulting receiver-side wavefield against the
%    saved source-side wavefield (from the matching forward run) every T
%    time steps to accumulate the RTM image contribution for one shot.
%
%  Inputs:
%    dx, dz   - Grid spacing in x and z [m]
%    nt       - Number of time samples
%    dt       - Time step [s]
%    v        - Background (smoothed) velocity model [nz x nx], PML-padded
%    xr, zr   - Receiver grid indices for this shot
%    data     - Scattered-data residual for this shot
%    d_pml    - PML thickness in grid cells
%    ss       - Shot index (used to load the matching saved forward
%               wavefield from ./Wavefields/UU_<ss>.mat)
%    Auu...Aww- MNAD operator matrices for the BACKGROUND velocity model
%               (same as those used to generate the saved
%               forward wavefield UU for this shot)
%
%  Output:
%    Im       - Migration image contribution for this shot (vectorized,
%               interior region only, PML stripped), to be summed/stacked
%               across shots 
% ========================================================================

% % --------- Load the saved source-side wavefield snapshots for this shot
load(['./Wavefields/UU_' num2str(ss)])

[nz,nx] = size(v) ;
V2dt2 = (v.*dt).^2 ;

% % --------- Snapshot storage sizing (must match the forward run's T/nf) 
% Recomputed here from the same background model v that produced UU, so
% this should reproduce the same snapshot interval T used when UU was saved
nf = floor(nt/(2*3*max(v(:))/min(v(:)))) ;
T = floor(nt/nf/2) ;
ttt = 1 ;
Im = 0 ;

% % --------- NOTE: assumes `data` is a cell {datau, datav, dataw} as returned by
% fwd_PDE_exact_MNAD, Here we used data{1} = datau but you may using
% data{2} = datav and  data{3} = dataw while the source-side fields also 
% must be horizontal or vertical gradient field, respectively. 
datau = data{1} ;


% % --------- Reconstructed receiver-side wavefield state 
pU_r = zeros(nz,nx) ;       % previous wavefield in time: U^(t-1)
U_r = zeros(nz,nx) ;        % current wavefield in time: U^(t)

pV = zeros(nz,nx) ;         % previous horizontal gradient field in time: V^(t-1)
V = zeros(nz,nx) ;          % current horizontal gradient field in time: V^(t)

pW = zeros(nz,nx) ;         % previous vertical gradient field in time: W^(t-1)
W = zeros(nz,nx) ;          % current vertical gradient field in time: W^(t)


% % --------- PML auxiliary variables 
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
R = 10^-4 ;
[c1xl,c2xl,c3xl,c1xr,c2xr,c3xr,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb] = C_cof(R,d_pml,dx,dz,dt,nx,nz,v) ;

% % --------- Time-reversed propagation loop
for tt = nt:-1:1
    
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
           
    % --------- Advance the reconstructed receiver-side wavefield one step backward
    nU_r = reshape(Auu*U_r(:),nz,nx) + reshape(Auv*V(:),nz,nx) + reshape(Auw*W(:),nz,nx) - pU_r - V2dt2.*(U1_x + U1_z) ;
    nV = reshape(Avv*V(:),nz,nx) + reshape(Avu*U_r(:),nz,nx) + reshape(Avw*W(:),nz,nx) - pV - V2dt2.*(V1_x + V1_z) ;
    nW = reshape(Aww*W(:),nz,nx) + reshape(Awu*U_r(:),nz,nx) + reshape(Awv*V(:),nz,nx) - pW - V2dt2.*(W1_x + W1_z) ;

    % --------- Inject scattered-data residual at receiver locations
    for rr = 1:length(xr)
        U_r(zr(rr),xr(rr)) = U_r(zr(rr),xr(rr)) + V2dt2(zr(rr),xr(rr)).*datau(tt,rr) ;
    end


    %     if mod(tt,50)==0, imagesc(U); drawnow; end
    
    pU_r = U_r ;
    U_r = nU_r ;
    
    pV = V ;
    V = nV ;
    
    pW = W ;
    W = nW ;
    
    
    % --------- cross-correlate the wavefields every T steps (imaging condition)
    if rem(tt,T) == 0
        US = UU(:,end-ttt+1) ;                                % source-side snapshot, reverse-indexed
        internal = U_r(d_pml+1:end-d_pml,d_pml+1:end-d_pml) ; % receiver-side field, interior only
        UR = internal(:) ;
        
        Im = Im + US.*UR ;
        ttt = ttt + 1 ;
    end

end

Im = T*Im ;

end
