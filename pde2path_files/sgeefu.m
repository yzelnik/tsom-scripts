function [p,idx]=sgeefu(p,u)
E=abs(p.mat.K*u(1:p.np)); 
[po,tri,ed]=getpte(p); E=pdeintrp(po,tri,E);
p.sol.err=max(max(E)); idx=p.pdeo.selectElements2Refine(E,p.nc.sig); 

