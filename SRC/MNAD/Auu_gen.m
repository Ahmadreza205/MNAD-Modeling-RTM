function [A] = Auu_gen(dx,dz,dt,v)

% % ========================================================================
%  MNAD Stencil Operator Generator: Auu (u <- u coupling)
% ========================================================================
%  Purpose:
%    Builds the sparse operator matrix Auu used in the
%    MNAD equation for wavefield propagation.
%    Auu uses 9-point stencil (4 horizontal neighbors +
%    4 vertical neighbors + self term). Coefficients follow Appendix A
%    of the manuscript.
%
%  Inputs:
%    dx, dz - Grid spacing in x and z [m]
%    dt     - Time step [s]
%    v      - Velocity model [nz x nx] (PML-padded)
%
%  Output:
%    A      - Sparse [nx*nz x nx*nz] stencil operator matrix Auu
% ========================================================================
[nz,nx] = size(v) ;
idx = reshape(1:nx*nz,nz,nx) ;    % Linear index lookup: idx(j,i) -> grid point (z=j, x=i)

% % u(x,z)
% u(i,j+1)
j = 2:nz ;
i = 1:nx ;
I1 = idx(j,i) ;
J1 = idx(j-1,i) ;
v1 = (24/11)*v(j,i).^2*dt^2/dz^2 - v(j,i).^4*dt^4*(15/(11*dz^4) + 1/(dx^2*dz^2))  ;

% u(i+1,j)
j = 1:nz ;
i = 2:nx ;
I2 = idx(j,i) ;
J2 = idx(j,i-1) ;
v2 = (24/11)*v(j,i).^2*dt^2/dx^2 - v(j,i).^4*dt^4*(15/(11*dx^4) + 1/(dx^2*dz^2)) ;

% u(i,j)
j = 1:nz ;
i = 1:nx ;
I3 = idx(j,i) ;
J3 = idx(j,i) ;
v3 = 2 - (48/11)*v(j,i).^2*dt^2*(1/dx^2 + 1/dz^2) + (30/11)*v(j,i).^4*dt^4*(1/dx^4 + 1/dz^4 + 11/(15*dx^2*dz^2)) ;

% u(i-1,j)
j = 1:nz ;
i = 1:nx-1 ;
I4 = idx(j,i) ;
J4 = idx(j,i+1) ;
v4 = (24/11)*v(j,i).^2*dt^2/dx^2 - v(j,i).^4*dt^4*(15/(11*dx^4) + 1/(dx^2*dz^2)) ;

% u(i,j-1)
j = 1:nz-1 ;
i = 1:nx ;
I5 = idx(j,i) ;
J5 = idx(j+1,i) ;
v5 = (24/11)*v(j,i).^2*dt^2/dz^2 - v(j,i).^4*dt^4*(15/(11*dz^4) + 1/(dx^2*dz^2)) ;

% u(i+1,j+1)
j = 2:nz ;
i = 2:nx ;
I6 = idx(j,i) ;
J6 = idx(j-1,i-1) ;
v6 = v(j,i).^4*dt^4/(2*dx^2*dz^2) ;

% u(i+1,j-1)
j = 1:nz-1 ;
i = 2:nx ;
I7 = idx(j,i) ;
J7 = idx(j+1,i-1) ;
v7 = v(j,i).^4*dt^4/(2*dx^2*dz^2) ;

% u(i-1,j+1)
j = 2:nz ;
i = 1:nx-1 ;
I8 = idx(j,i) ;
J8 = idx(j-1,i+1) ;
v8 = v(j,i).^4*dt^4/(2*dx^2*dz^2) ;

% u(i-1,j-1)
j = 1:nz-1 ;
i = 1:nx-1 ;
I9 = idx(j,i) ;
J9 = idx(j+1,i+1) ;
v9 = v(j,i).^4*dt^4/(2*dx^2*dz^2) ;

A = sparse([I1(:);I2(:);I3(:);I4(:);I5(:);I6(:);I7(:);I8(:);I9(:)],...
[J1(:);J2(:);J3(:);J4(:);J5(:);J6(:);J7(:);J8(:);J9(:)],...
[v1(:);v2(:);v3(:);v4(:);v5(:);v6(:);v7(:);v8(:);v9(:)],...
nx*nz,nx*nz) ;
end
