function [U1_z,U1_zt,U2_zt,U3_zt,U1_zb,U2_zb,U3_zb] =...
    Uz123_gen(U1_zt,U2_zt,U3_zt,U1_zb,U2_zb,U3_zb,U,dz,nz,nx,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb)

U_az = zeros(nz,nx) ;
U_az(1:nz-1,:) = U(2:nz,:) ;
U_dz = zeros(nz,nx) ;
U_dz(2:nz,:) = U(1:nz-1,:) ;

U3_zt = c1zt .* U3_zt + (c2zt./(2*dz)) .* (U_dz - U_az) ;
U2_zt = c1zt .* U2_zt + (c2zt./dz) .* (U_dz - U_az) + c3zt .* ((U_dz - 2.*U + U_az)./(dz.^2) - U3_zt) ;
U1_zt = c1zt .* U1_zt + (c2zt./(2*dz)) .* (U_dz - U_az) + c3zt .* (2.*(U_dz - 2.*U + U_az)./(dz.^2) - U2_zt) ;

U3_zb = c1zb .* U3_zb + (c2zb./(2*dz)) .* (U_az - U_dz) ;
U2_zb = c1zb .* U2_zb + (c2zb./dz) .* (U_az - U_dz) + c3zb .* ((U_az - 2.*U + U_dz)./(dz.^2) - U3_zb) ;
U1_zb = c1zb .* U1_zb + (c2zb./(2*dz)) .* (U_az - U_dz) + c3zb .* (2.*(U_az - 2.*U + U_dz)./(dz.^2) - U2_zb) ;

U1_z = U1_zb + U1_zt ;
end