function [c1xl,c2xl,c3xl,c1xr,c2xr,c3xr,c1zt,c2zt,c3zt,c1zb,c2zb,c3zb] = C_cof(R,d_pml,dx,dz,dt,nx,nz,v)

%%%%%%%%%%%%%% X accesories %%%%%%%%%%%%%%
D_pml = d_pml *dx ;

x_b = [D_pml:-dx:dx,zeros(1,nx-2*d_pml),dx:dx:D_pml] ; 
X = repmat(x_b,nz,1) ;

% Left
V_xl = v;
V_xl(:,d_pml+1:nx) = 0 ;

sig_xl = ((3*log(1/R))/(2*D_pml^3)).*V_xl ;
sigxl = sig_xl.*(X.^2) ;
dsigxl = 2.*sig_xl .* X ;
denxl = (1+sigxl.*(dt/2)) ;

c1xl = sparse((1-sigxl.*(dt/2)) ./ denxl);
c2xl = sparse((dsigxl .*dt) ./ denxl) ;
c3xl = sparse((sigxl .*dt) ./ denxl) ;

% Right
V_xr = v;
V_xr(:,1:nx-d_pml) = 0 ;

sig_xr = ((3*log(1/R))/(2*D_pml^3)).*V_xr ;
sigxr = sig_xr.*(X.^2) ;
dsigxr = 2.*sig_xr .* X ;
denxr = (1+sigxr.*(dt/2)) ;

c1xr = sparse((1-sigxr.*(dt/2)) ./ denxr) ;
c2xr = sparse((dsigxr .*dt) ./ denxr) ;
c3xr = sparse((sigxr .*dt) ./ denxr) ;

%%%%%%%%%%%%%% Z accesories %%%%%%%%%%%%%%
D_pml = d_pml *dz ;

z_b = [D_pml:-dz:dz,zeros(1,nz-2*d_pml),dz:dz:D_pml] ; 
Z = repmat(z_b',1,nx) ;

% Top
V_zt = v;
V_zt(d_pml+1:nz,:) = 0 ;

sig_zt = ((3*log(1/R))/(2*D_pml^3)).*V_zt ;
sigzt = sig_zt.*(Z.^2) ;
dsigzt = 2.*sig_zt .* Z ;
denzt = (1+sigzt.*(dt/2)) ;

c1zt = sparse((1-sigzt.*(dt/2)) ./ denzt) ;
c2zt = sparse((dsigzt .*dt) ./ denzt) ;
c3zt = sparse((sigzt .*dt) ./ denzt) ;

% Bottom
V_zb = v;
V_zb(1:nz-d_pml,:) = 0 ;

sig_zb = ((3*log(1/R))/(2*D_pml^3)).*V_zb ;
sigzb = sig_zb.*(Z.^2) ;
dsigzb = 2.*sig_zb .* Z ;
denzb = (1+sigzb.*(dt/2)) ;

c1zb = sparse((1-sigzb.*(dt/2)) ./ denzb) ;
c2zb = sparse((dsigzb .*dt) ./ denzb) ;
c3zb = sparse((sigzb .*dt) ./ denzb) ;

end
