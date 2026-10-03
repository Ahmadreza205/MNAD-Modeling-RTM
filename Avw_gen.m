function [A] = Avw_gen(dx,dz,dt,v)

% % ========================================================================
%  MNAD Stencil Operator Generator: Avw (w <- v coupling)
% ========================================================================
%  Purpose:
%    Builds the sparse operator matrix Avw used in the
%    MNAD equation for horizontal gradient field propagation.
%    Avw uses 9-point stencil (4 horizontal neighbors +
%    4 vertical neighbors + self term). Coefficients follow Appendix A
%    of the manuscript.
%
%  Inputs:
%    dx, dz - Grid spacing in x and z [m]
%    dt     - Time step [s]
%    v      - Velocity model [nz x nx] (PML-padded)
%
%  Output:
%    A      - Sparse [nx*nz x nx*nz] stencil operator matrix Avw
% ========================================================================
[nz,nx] = size(v) ;
idx = reshape(1:nx*nz,nz,nx) ;    % Linear index lookup: idx(j,i) -> grid point (z=j, x=i)

% w(i+1,j+1)
j = 2:nz ;
i = 2:nx ;
I1 = idx(j,i) ;
J1 = idx(j-1,i-1) ;
v1 = v(j,i).^4*dt^4/(4*dz^3*dx) - v(j,i).^2*dt^2/(4*dz*dx) ;

% w(i+1,j-1)
j = 1:nz-1 ;
i = 2:nx ;
I2 = idx(j,i) ;
J2 = idx(j+1,i-1) ;
v2 = -v(j,i).^4*dt^4/(4*dz^3*dx) + v(j,i).^2*dt^2/(4*dz*dx) ;

% w(i-1,j+1)
j = 2:nz ;
i = 1:nx-1 ;
I3 = idx(j,i) ;
J3 = idx(j-1,i+1) ;
v3 = -v(j,i).^4*dt^4/(4*dz^3*dx) + v(j,i).^2*dt^2/(4*dz*dx) ;

% w(i-1,j-1)
j = 1:nz-1 ;
i = 1:nx-1 ;
I4 = idx(j,i) ;
J4 = idx(j+1,i+1) ;
v4 = v(j,i).^4*dt^4/(4*dz^3*dx) - v(j,i).^2*dt^2/(4*dz*dx) ;


A = sparse([I1(:);I2(:);I3(:);I4(:)],...
[J1(:);J2(:);J3(:);J4(:)],...
[v1(:);v2(:);v3(:);v4(:)],...
nx*nz,nx*nz) ;
end