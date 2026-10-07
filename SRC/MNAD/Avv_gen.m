function [A] = Avv_gen(dx,dz,dt,v)

% % ========================================================================
%  MNAD Stencil Operator Generator: Avv (v <- v coupling)
% ========================================================================
%  Purpose:
%    Builds the sparse operator matrix Avv used in the
%    MNAD equation for horizontal gradient field propagation.
%    Avv uses 9-point stencil (4 horizontal neighbors +
%    4 vertical neighbors + self term). Coefficients follow Appendix A
%    of the manuscript.
%
%  Inputs:
%    dx, dz - Grid spacing in x and z [m]
%    dt     - Time step [s]
%    v      - Velocity model [nz x nx] (PML-padded)
%
%  Output:
%    A      - Sparse [nx*nz x nx*nz] stencil operator matrix Avv
% ========================================================================
[nz,nx] = size(v) ;
idx = reshape(1:nx*nz,nz,nx) ;    % Linear index lookup: idx(j,i) -> grid point (z=j, x=i)

% v(i,j+1)
j = 2:nz ;
i = 1:nx ;
I1 = idx(j,i) ;
J1 = idx(j-1,i) ;
v1 = v(j,i).^2*dt^2/dz^2 - 2*v(j,i).^4*dt^4/(dx^2*dz^2) ;

% v(i+1,j)
j = 1:nz ;
i = 2:nx ;
I2 = idx(j,i) ;
J2 = idx(j,i-1) ;
v2 = -(3/2)*v(j,i).^2*dt^2/dx^2 + v(j,i).^4*dt^4*(5/(2*dx^4) + 1/(2*dx^2*dz^2)) ;

% v(i,j)
j = 1:nz ;
i = 1:nx ;
I3 = idx(j,i) ;
J3 = idx(j,i) ;
v3 = 2 - v(j,i).^2*dt^2*(12/dx^2 + 2/dz^2) + v(j,i).^4*dt^4*(10/dx^4 + 4/(dx^2*dz^2)) ;

% v(i-1,j)
j = 1:nz ;
i = 1:nx-1 ;
I4 = idx(j,i) ;
J4 = idx(j,i+1) ;
v4 = -(3/2)*v(j,i).^2*dt^2/dx^2 + v(j,i).^4*dt^4*(5/(2*dx^4) + 1/(2*dx^2*dz^2)) ;

% v(i,j-1)
j = 1:nz-1 ;
i = 1:nx ;
I5 = idx(j,i) ;
J5 = idx(j+1,i) ;
v5 = v(j,i).^2*dt^2/dz^2 - 2*v(j,i).^4*dt^4/(dx^2*dz^2) ;

% v(i+1,j+1)
j = 2:nz ;
i = 2:nx ;
I6 = idx(j,i) ;
J6 = idx(j-1,i-1) ;
v6 = -v(j,i).^4*dt^4/(4*dx^2*dz^2) ;

% v(i+1,j-1)
j = 1:nz-1 ;
i = 2:nx ;
I7 = idx(j,i) ;
J7 = idx(j+1,i-1) ;
v7 = -v(j,i).^4*dt^4/(4*dx^2*dz^2) ;

% v(i-1,j+1)
j = 2:nz ;
i = 1:nx-1 ;
I8 = idx(j,i) ;
J8 = idx(j-1,i+1) ;
v8 = -v(j,i).^4*dt^4/(4*dx^2*dz^2) ;

% v(i-1,j-1)
j = 1:nz-1 ;
i = 1:nx-1 ;
I9 = idx(j,i) ;
J9 = idx(j+1,i+1) ;
v9 = -v(j,i).^4*dt^4/(4*dx^2*dz^2) ;

A = sparse([I1(:);I2(:);I3(:);I4(:);I5(:);I6(:);I7(:);I8(:);I9(:)],...
[J1(:);J2(:);J3(:);J4(:);J5(:);J6(:);J7(:);J8(:);J9(:)],...
[v1(:);v2(:);v3(:);v4(:);v5(:);v6(:);v7(:);v8(:);v9(:)],...
nx*nz,nx*nz) ;
end
