function [A] = Auv_gen(dx,dz,dt,v)
% % ========================================================================
%  MNAD Stencil Operator Generator: Auv (v <- u coupling)
% ========================================================================
%  Purpose:
%    Builds the sparse operator matrix Auv used in the
%    MNAD equation for wavefield propagation.
%    Auv uses 9-point stencil (4 horizontal neighbors +
%    4 vertical neighbors + self term). Coefficients follow Appendix A
%    of the manuscript.
%
%  Inputs:
%    dx, dz - Grid spacing in x and z [m]
%    dt     - Time step [s]
%    v      - Velocity model [nz x nx] (PML-padded)
%
%  Output:
%    A      - Sparse [nx*nz x nx*nz] stencil operator matrix Auv
% ========================================================================
[nz,nx] = size(v) ;
idx = reshape(1:nx*nz,nz,nx) ;    % Linear index lookup: idx(j,i) -> grid point (z=j, x=i)


% v(i+1,j)
j = 1:nz ;
i = 2:nx ;
I1 = idx(j,i) ;
J1 = idx(j,i-1) ;
v1 = (-20/33)*v(j,i).^2*dt^2/dx + v(j,i).^4*dt^4*(47/(66*dx^3) +1/(6*dx*dz^2)) ;

% v(i-1,j)
j = 1:nz ;
i = 1:nx-1 ;
I2 = idx(j,i) ;
J2 = idx(j,i+1) ;
v2 = (20/33)*v(j,i).^2*dt^2/dx - v(j,i).^4*dt^4*(47/(66*dx^3) +1/(6*dx*dz^2)) ;

% v(i+1,j+1)
j = 2:nz ;
i = 2:nx ;
I3 = idx(j,i) ;
J3 = idx(j-1,i-1) ;
v3 = -v(j,i).^4*dt^4/(12*dx*dz^2) ;

% v(i+1,j-1)
j = 1:nz-1 ;
i = 2:nx ;
I4 = idx(j,i) ;
J4 = idx(j+1,i-1) ;
v4 = -v(j,i).^4*dt^4/(12*dx*dz^2) ;

% v(i-1,j+1)
j = 2:nz ;
i = 1:nx-1 ;
I5 = idx(j,i) ;
J5 = idx(j-1,i+1) ;
v5 = v(j,i).^4*dt^4/(12*dx*dz^2) ;

% v(i-1,j-1)
j = 1:nz-1 ;
i = 1:nx-1 ;
I6 = idx(j,i) ;
J6 = idx(j+1,i+1) ;
v6 = v(j,i).^4*dt^4/(12*dx*dz^2) ;

% v(i+2,j)
j = 1:nz ;
i = 3:nx ;
I7 = idx(j,i) ;
J7 = idx(j,i-2) ;
v7 = (1/132)*v(j,i).^2*dt^2/dx - (1/66)*v(j,i).^4*dt^4/dx^3 ;

% v(i-2,j)
j = 1:nz ;
i = 1:nx-2 ;
I8 = idx(j,i) ;
J8 = idx(j,i+2) ;
v8 = (-1/132)*v(j,i).^2*dt^2/dx + (1/66)*v(j,i).^4*dt^4/dx^3 ;

A = sparse([I1(:);I2(:);I3(:);I4(:);I5(:);I6(:);I7(:);I8(:)],...
[J1(:);J2(:);J3(:);J4(:);J5(:);J6(:);J7(:);J8(:)],...
[v1(:);v2(:);v3(:);v4(:);v5(:);v6(:);v7(:);v8(:)],...
nx*nz,nx*nz) ;
end