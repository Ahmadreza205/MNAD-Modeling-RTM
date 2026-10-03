function [U1_x,U1_xl,U2_xl,U3_xl,U1_xr,U2_xr,U3_xr] =...
    Ux123_gen(U1_xl,U2_xl,U3_xl,U1_xr,U2_xr,U3_xr,U,dx,nz,nx,c1xl,c2xl,c3xl,c1xr,c2xr,c3xr)

U_ax = zeros(nz,nx) ;
U_ax(:,1:nx-1) = U(:,2:nx) ;
U_dx = zeros(nz,nx) ;
U_dx(:,2:nx) = U(:,1:nx-1) ;

U3_xl = c1xl .* U3_xl + (c2xl./(2*dx)) .* (U_dx - U_ax) ;
U2_xl = c1xl .* U2_xl + (c2xl./dx) .* (U_dx - U_ax) + c3xl .* ((U_dx - 2.*U + U_ax)./(dx.^2) - U3_xl) ;
U1_xl = c1xl .* U1_xl + (c2xl./(2*dx)) .* (U_dx - U_ax) + c3xl .* (2.*(U_dx - 2.*U + U_ax)./(dx.^2) - U2_xl) ;

U3_xr = c1xr .* U3_xr + (c2xr./(2*dx)) .* (U_ax - U_dx) ;
U2_xr = c1xr .* U2_xr + (c2xr./dx) .* (U_ax - U_dx) + c3xr .* ((U_ax - 2.*U + U_dx)./(dx.^2) - U3_xr) ;
U1_xr = c1xr .* U1_xr + (c2xr./(2*dx)) .* (U_ax - U_dx) + c3xr .* (2.*(U_ax - 2.*U + U_dx)./(dx.^2) - U2_xr) ;

U1_x = U1_xl + U1_xr  ;

end