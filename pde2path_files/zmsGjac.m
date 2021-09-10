function Gu=zmsGjac(p,u)
par=u(p.nu+1:end); n=p.np; 
q=par(1); nu=par(2); al=par(3); eta=par(4); rp=par(5); f=par(6); 
dh=par(7); dw=par(8); pp=par(9);  gam=par(10); 
%spaf=exp(-p.pdeo.grid.p(1,:)'.^2/50)*f/1000+f;
if(isempty(who('kf'))) kf=0.15343097; end;
spaf=f*(1+0.5*(1+cos((kf)*p.pdeo.grid.p(1,:)')));
%spaf=f*(1+0.5*(1+cos((1.1*0.43)*p.pdeo.grid.p(1,:)')));

%u=splitu(p,u); b=u(:,1); w=u(:,2); h=u(:,3); 
b=u(1:p.np); w=u(p.np+1:2*p.np); h=u(2*p.np+1:3*p.np); 
gb=nu*w.*(1+eta*b).^2; gw=gam*b.*(1+eta*b).^2; 
I=al*(b+q*spaf)./(b+q); ov=ones(n,1); zv=0*ov; 
gbb=2*nu*eta*w.*(1+eta*b); gbw=nu*(1+eta*b).^2; 
gwb=gam*(1+eta*b).^2+2*gam*eta*b.*(1+eta*b); 
Ib=al*q*(1-spaf)./((b+q).^2); 
f1b=gbb.*b.*(1-b)-1+gb.*(1-2*b); f1w=gbw.*b.*(1-b); f1h=zv; 
f2b=Ib.*h+nu*w*rp-gwb.*w; f2w=-nu*(1-rp*b)-gw; f2h=I; 
f3b=-Ib.*h; f3w=zv;  f3h=-I; 
K=p.mat.K; M=p.mat.M(1:p.np,1:p.np); 
hK=K*spdiags(h,0,n,n); % hK=bsxfun(@times,K',h)'; 
Gu=[[K-M*spdiags(f1b,0,n,n),-M*spdiags(f1w,0,n,n),-M*spdiags(f1h,0,n,n) ];
   [-M*spdiags(f2b,0,n,n),dw*K-M*spdiags(f2w,0,n,n), -M*spdiags(f2h,0,n,n)]
   [-M*spdiags(f3b,0,n,n),-M*spdiags(f3w,0,n,n), 2*dh*hK-M*spdiags(f3h,0,n,n)]]; 
end      
