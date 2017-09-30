function result = calcV( k,x1,x2,a,b,cp,strike )
% Computes the V_k coefficients
%   Detailed explanation goes here

[chi,psi] = coeff_b(k,x1,x2,a,b);
result = 2.*cp.*strike./(b-a).*(chi-psi);

%NOTE: cp = -1 for put

end

