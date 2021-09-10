function VsOut=L_GrzMng(Vs,Ps,Es)
% Grazing managment model - Local terms

if(~isfield(Es,'JacMode'))
   Es.JacMode=0;
end;

% Initialization
B=Vs(:,1);
W=Vs(:,2);
if(Es.JacMode==0)      % Model equations

    dB = Ps.lambda.*B.*W.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa) - Ps.mu.*B;
    dW = Ps.P - Ps.nu.*W./(1 + Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*W.*(1 + Ps.eta.*B).^2;
    VsOut = [dB,dW];
else               % Jacobian of equations
    BdB = (-Ps.kappa.* Ps.mu + (1 + B.* Ps.eta).* (Ps.kappa + B.* (-2 - 4* B.* Ps.eta + 3* Ps.eta.* Ps.kappa)).* Ps.lambda.* W)./Ps.kappa;
    BdW = Ps.lambda.*B.*(1 + Ps.eta.*B).^2.*(1 - B./Ps.kappa);
    WdB = - Ps.gamma.*W.*(1 + 4* B.* Ps.eta + 3* B.^2.* Ps.eta.^2) - Ps.rho.*Ps.nu.*W./Ps.kappa .* ((1 + Ps.rho.*B./Ps.kappa).^(-2));
    WdW = - Ps.nu./(1 + Ps.rho.*B./Ps.kappa) - Ps.gamma.*B.*(1 + Ps.eta.*B).^2 ;

    % written in a large sparse matrix format
    VsOut = ArrangeJacobian([BdB BdW;WdB WdW],Ps,Es);
end;

end
