
function [c1, c2, c4] = getCumulants(T,r,d,model)

if strcmp(model.ID,'BlackScholes')
    BS_sigma = model.params(1);
    c1 = (r-d-0.5*BS_sigma^2)*T;
    c2 = BS_sigma^2*T;
    c4 = 0.0;
    
    
elseif strcmp(model.ID,'Merton')
    BS_sigma = sqrt( model.params(1)^2 + model.params(3)^2* model.params(4)  );
    c1 = (r-d-0.5*BS_sigma^2)*T;
    c2 = BS_sigma^2*T;
    c4 = 0.0;
    
elseif strcmp(model.ID,'Kou')
    Kou_sigma = sqrt( model.params(1) * model.params(1) + 2*model.params(2) *( model.params(3)./ (model.params(4))^2 + ( 1-model.params(3))./ (model.params(5))^2  )  );
    c1 = (r-d-0.5*Kou_sigma^2 )*T ;
    
    delta  =  model.params(2) *(     model.params(3)./  model.params(4)  + ( 1-model.params(3))./ model.params(5)  ...
              - (   model.params(3) .* model.params(4) ./  ( model.params(4)-1)  + ( 1-model.params(3)) .* model.params(5) ./ ( 1+ model.params(5) )  - 1 )   );    
    c1 = ( r-d-0.5*model.params(1) * model.params(1) + delta )*T ; 
    
    c2 = Kou_sigma^2*T;
    c4 = 0.0;
    
    gamma1 = r- d- 1/2*model.params(1) * model.params(1) - model.params(2)*(  model.params(3)./ (model.params(4)-1) + ( 1-model.params(3))./ (model.params(5)+1)   );
    %% Edited by Lian 2016 Jun in Berkely
    c1 = T*( gamma1 + model.params(2) *( model.params(3)./  model.params(4) + ( 1-model.params(3))./ model.params(5) ) );
    c2 = T*(  model.params(1) * model.params(1) + 2*model.params(2)*(  model.params(3)./  model.params(4)^2 + ( 1-model.params(3))./ model.params(5)^2  ) );
    c4 = 24*T*model.params(2)*(  model.params(3)./  model.params(4)^4 + ( 1-model.params(3))./ model.params(5)^4  ) ;
    
elseif strcmp(model.ID,'Heston')
    vInst = model.params(1); vLong = model.params(2); kappa = model.params(3);
    omega = model.params(4); rho = model.params(5);
    c1 = (r-d)*T + 0.5*((1-exp(-kappa*T))*(vLong-vInst)/kappa - vLong*T);
    c2 = 1/8/kappa^2*(omega*T*exp(-kappa*T)*(vInst-vLong)*(8*kappa*rho-4*omega)...
        +8*rho*omega*(1-exp(-kappa*T))*(2*vLong-vInst)...
        +2*vLong*T*(-4*kappa*rho*omega + omega^2 + 4*kappa^2)...
        +omega^2/kappa*((vLong-2*vInst)*exp(-2*kappa*T)+vLong*(6*exp(-kappa*T)-7)+2*vInst)...
        +8*kappa*(vInst-vLong)*(1-exp(-kappa*T)));
    c4 = 0.0;
elseif strcmp(model.ID,'VarianceGamma')
    VG_sigma = model.params(1); VG_nu = model.params(2); VG_theta = model.params(3);
    omega = 1/VG_nu * log(1 - VG_theta * VG_nu - 0.5 * VG_nu * VG_sigma^2);
    
    c1 = (r-d+omega+VG_theta)*T;
    c2 = (VG_sigma^2+VG_nu*VG_theta^2)*T;
    c4 = 3*(VG_sigma^4*VG_nu+2*VG_theta^4*VG_nu^3+4*VG_sigma^2*VG_theta^2*VG_nu^2)*T;
elseif strcmp(model.ID,'CGMY')
    CGMY_C = model.params(1); CGMY_G = model.params(2); CGMY_M = model.params(3); CGMY_Y = model.params(4);
    c1 = (r-d)*T + CGMY_C*T* gamma(1-CGMY_Y)*(CGMY_M^(CGMY_Y-1)-CGMY_G^(CGMY_Y-1));
    c2 = CGMY_C*T*gamma(2-CGMY_Y)*(CGMY_M^(CGMY_Y-2)+CGMY_G^(CGMY_Y-2));
    c4 = CGMY_C*T*gamma(4-CGMY_Y)*(CGMY_M^(CGMY_Y-4)+CGMY_G^(CGMY_Y-4));
elseif strcmp(model.ID,'NIG')
    alphaNIG = model.params(1); betaNIG = model.params(2); deltaNIG = model.params(3);
    mu = r - d;
    c1 = ((betaNIG*deltaNIG)/(sqrt(alphaNIG^2 - betaNIG^2)) + mu) * T;
    c2 = T * alphaNIG^2 * deltaNIG / ((alphaNIG^2 - betaNIG^2)^(3/2));
    c4 = 3 * T * alphaNIG^2 * (alphaNIG^2 + 4*betaNIG^2) * deltaNIG / ((alphaNIG^2 - betaNIG^2)^(7/2));
elseif strcmp(model.ID,'Bates')
    vInst = model.params(1); vLong = model.params(2); kappa = model.params(3);
    omega = model.params(4); rho = model.params(5); lambda = model.params(6);
    muj = model.params(7); sigmaj = model.params(8);
    c1 = (r-d).*T + 0.5*((1-exp(-kappa*T))*(vLong-vInst)/kappa - vLong*T) + lambda*muj*T;
    c2 = 1/8/kappa^2*(omega*T.*exp(-kappa*T)*(vInst-vLong)*(8*kappa*rho-4*omega)...
        +8*rho*omega*(1-exp(-kappa*T))*(2*vLong-vInst)...
        +2*vLong*T*(-4*kappa*rho*omega + omega^2 + 4*kappa^2)...
        +omega^2/kappa*((vLong-2*vInst)*exp(-2*kappa*T)...
        +vLong*(6*exp(-kappa*T)-7)+2*vInst)...
        +8*kappa*(vInst-vLong)*(1-exp(-kappa*T))) + lambda*(sigmaj^2+muj^2)*T;
    c4 = 0.0 + lambda * (3*sigmaj^2*(sigmaj^2+2*muj^2)+muj^4).*T;
    
end

end
