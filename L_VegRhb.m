function VsOut=L_VegRhb(Vs,Ps,Es)
% Vegetation rehabilitation model - Local terms

% Initialization
B=Vs(:,1); 
W=Vs(:,2); 
H=Vs(:,3); 

if(Es.JacMode==0)      % Model equations

    dB = Ps.lambda.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu*B;
    dW = Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H - Ps.nu*W.*(1 - Ps.rho.*B./Ps.kappa) - Ps.gamma.*W.*(1 + Ps.eta.*B).^2.*B;
    dH = Ps.P - Ps.alpha.*(B + Ps.q.*Ps.f)./(B + Ps.q).*H;

    VsOut = [dB,dW,dH];
else                % Jacobian of equations not implemented
    
    % written in a large sparse matrix format 
    VsOut = 0;%ArrangeJacobian([BdB BdW BdH;WdB WdW WdH; HdB HdW HdH],Ps,Es);
end;


