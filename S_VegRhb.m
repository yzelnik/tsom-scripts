function Out=S_VegRhb(Vs,Ps,Es)
% Vegetation rehabilitation model - Spatial function

if(~isfield(Ps,'Nld'))
    Ps.Nld=0;
end;
if(~isfield(Ps,'H2') || Ps.H2==0)
   Ps.H2=Ps.VarNum; % By default the H^2 field is assumed to be the last one
end;
if(size(Ps.Ds,2)==1) Ps.Ds=Ps.Ds'; end; % Make sure Ds is in a row shape

if(isfield(Es,'SetupMode') && Es.SetupMode)
   % Pre calculate spatial matrix, for future use
   Out = Ps;
   Out.Derv2Mat = DervSM(2,Ps,Es);
   if(Ps.Nld==0) % More direct way, just apply second derivative on H^2

   else          % Less direct, calculate SM including values of H
       Es.JacMode=1;
       Es.SetupMode=0;
       Out.SpaMat=S_VegRhb(Vs,Out,Es);
       Es.JacMode=0;    
   end;
else            % Normal run
   len=Ps.Nx*Ps.Ny; 
   
   if(~isfield(Ps,'Derv2Mat'))        % Caclculate spatial matrix if needed
        Ps.Derv2Mat = DervSM(2,Ps,Es);
   end;
   
   if(~isfield(Es,'JacMode') || (Es.JacMode==0))	% Model equations
        if(Ps.Nld==0)
            Out = [];
            for ind=1:Ps.VarNum
                if(ind==Ps.H2)
                    tmpout = Ps.Derv2Mat*(Vs(:,ind).^2).*Ps.Ds(:,ind);
                else
                    tmpout = Ps.Derv2Mat*Vs(:,ind).*Ps.Ds(:,ind);
                end;
            Out  = [Out tmpout];
            end;
        else % Use the jacobian as the SM
            Es.JacMode=1;
            Out=S_VegRhb(Vs,Ps,Es);
            Es.JacMode=0;
        end;
   else			% Jacobian of equations
       
        Out = sparse(len*Ps.VarNum,len*Ps.VarNum);
        for ind=1:Ps.VarNum
            if(ind==Ps.H2)
                H     = Vs(:,ind);
                grad  = GradSM(1,Ps,Es);
                part1 = sparse(H*ones(1,len)).*Ps.Derv2Mat;		% H0 * Derv2 (dH)
                part2 = sparse(diag(Ps.Derv2Mat*H));			% dH * Derv2 (H0)
                part3 = sparse(len,len);
                for jj=1:length(grad)  % 2* Derv1(H) * Derv1(dU)
                    part3=part3+2*sparse((grad{jj}*H)*ones(1,len)).*grad{jj};
                end;
                tmpout= 2* (part1 + part2 + part3);
            else
                tmpout= Ps.Derv2Mat;
            end;
            
            Out(len*(ind-1)+(1:len),len*(ind-1)+(1:len)) = tmpout*Ps.Ds(:,ind);
        end;
        
   end;
end;            % End normal run

end
