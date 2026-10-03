function [A] = Avu_gen(dx,dz,dt,v)

% % ========================================================================
%  MNAD Stencil Operator Generator: Avu (u <- v coupling)
% ========================================================================
%  Purpose:
%    Builds the sparse operator matrix Avu used in the
%    MNAD equation for horizontal gradient field propagation.
%    Avu uses 9-point stencil (4 horizontal neighbors +
%    4 vertical neighbors + self term). Coefficients follow Appendix A
%    of the manuscript.
%
%  Inputs:
%    dx, dz - Grid spacing in x and z [m]
%    dt     - Time step [s]
%    v      - Velocity model [nz x nx] (PML-padded)
%
%  Output:
%    A      - Sparse [nx*nz x nx*nz] stencil operator matrix Avu
% ========================================================================
[nz,nx] = size(v) ;
idx = reshape(1:nx*nz,nz,nx) ;    % Linear index lookup: idx(j,i) -> grid point (z=j, x=i)


% u(i+1,j)
j = 1:nz ;
i = 2:nx ;
I1 = idx(j,i) ;
J1 = idx(j,i-1) ;
v1 = -v(j,i).^4*dt^4*(15/(2*dx^5) - 1/(dx*dz^4) + 5/(2*dx^3*dz^2)) + v(j,i).^2*dt^2*(15/(2*dx^3) - 1/(dx*dz^2)) ;

% u(i-1,j)
j = 1:nz ;
i = 1:nx-1 ;
I2 = idx(j,i) ;
J2 = idx(j,i+1) ;
v2 = v(j,i).^4*dt^4*(15/(2*dx^5) - 1/(dx*dz^4) + 5/(2*dx^3*dz^2)) - v(j,i).^2*dt^2*(15/(2*dx^3) - 1/(dx*dz^2)) ;

% u(i+1,j+1)
j = 2:nz ;
i = 2:nx ;
I3 = idx(j,i) ;
J3 = idx(j-1,i-1) ;
v3 = v(j,i).^4*dt^4*(5/(4*dx^3*dz^2) - 1/(2*dx*dz^4)) + v(j,i).^2*dt^2/(2*dx*dz^2) ;

% u(i+1,j-1)
j = 1:nz-1 ;
i = 2:nx ;
I4 = idx(j,i) ;
J4 = idx(j+1,i-1) ;
v4 = v(j,i).^4*dt^4*(5/(4*dx^3*dz^2) - 1/(2*dx*dz^4)) + v(j,i).^2*dt^2/(2*dx*dz^2) ;

% u(i-1,j+1)
j = 2:nz ;
i = 1:nx-1 ;
I5 = idx(j,i) ;
J5 = idx(j-1,i+1) ;
v5 = -v(j,i).^4*dt^4*(5/(4*dx^3*dz^2) - 1/(2*dx*dz^4)) - v(j,i).^2*dt^2/(2*dx*dz^2) ;

% u(i-1,j-1)
j = 1:nz-1 ;
i = 1:nx-1 ;
I6 = idx(j,i) ;
J6 = idx(j+1,i+1) ;
v6 = -v(j,i).^4*dt^4*(5/(4*dx^3*dz^2) - 1/(2*dx*dz^4)) - v(j,i).^2*dt^2/(2*dx*dz^2) ;


A = sparse([I1(:);I2(:);I3(:);I4(:);I5(:);I6(:)],...
[J1(:);J2(:);J3(:);J4(:);J5(:);J6(:)],...
[v1(:);v2(:);v3(:);v4(:);v5(:);v6(:)],...
nx*nz,nx*nz) ;
end