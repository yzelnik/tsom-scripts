%% mod of tintx to explicit Euler 
%
%  [p,t1,ts,nc] = tintx(p,t0,ts,dt,nt,nc,pmod,smod)
function [p,t1,ts,nc]=expl_eu_tintx(p,t0,ts,dt,nt,nc,pmod,smod)
n=0; t=t0; [L,U,P,Q,R]=lu(p.mat.M);
while(n<nt) 
  G=pderesi(p,p.u); du=Q*(U\(L\(P*(R\G))));
  p.u(1:p.nu)=p.u(1:p.nu)-dt*du; 
  t=t+dt; n=n+1;
  if(mod(n,pmod)==0); 
      r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
      tits=['t=' mat2str(t,4) ', r=' mat2str(r,3)];
      plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      title(['u_1, ' tits],'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
      drawnow;
  end
  if(mod(n,smod)==0); 
    ps=p; p=[]; p.t=t;p.ts=ts; p.u=ps.u; 
    fname=[ps.file.pname,sprintf('%i',nc+n),'.mat']; save(fname,'p'); 
    p=ps;
  end
end 
t1=t; nc=nc+nt; 