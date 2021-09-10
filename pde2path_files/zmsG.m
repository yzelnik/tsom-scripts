function r=zmsG(p,u)
par=u(p.nu+1:end); 
q=par(1); nu=par(2); al=par(3); eta=par(4); rp=par(5); f=par(6); 
dh=par(7); dw=par(8); pp=par(9); gam=par(10); 
%spaf=exp(-p.pdeo.grid.p(1,:)'.^2/50)*f/1000+f;
if(isempty(who('kf'))) kf=0.15343097; end;
spaf=f*(1+0.5*(1+cos((kf)*p.pdeo.grid.p(1,:)')));
%spaf=f*(1+0.5*(1+cos((1.1*0.43)*p.pdeo.grid.p(1,:)')));

%u=splitu(p,u); 
b=u(1:p.np); w=u(p.np+1:2*p.np); h=u(2*p.np+1:3*p.np); 
gb=nu*w.*(1+eta*b).^2; gw=gam*b.*(1+eta*b).^2; 
I=al*(b+q*spaf)./(b+q);  
f1=gb.*b.*(1-b)-b; 
f2=I.*h-nu*w.*(1-rp.*b)-gw.*w; 
f3=pp-I.*h; 
K=p.mat.K; M=p.mat.M(1:p.np,1:p.np); 
r1=K*b-M*f1; r2=dw*K*w-M*f2; r3=dh*K*(h.^2)-M*f3; r=[r1;r2;r3]; 
end      
