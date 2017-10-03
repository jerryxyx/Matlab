function y = IV_test(r,T,S_0,cp,strike,Ngrid,Hes,c1,c2,c4)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

% Nstrike = size(strike,1);
% 
% x = repmat(double(log(S_0 ./ strike))',Ngrid,1);       % center

if T>=2
    L1 = 14;
    L2 = 14;
elseif T>=.1
%     L1 = 18;
%     L2 = 20;
       L1 = 12;
       L2 = 20;
else
%     L1 = 25;
%     L2 = 28;
       L1 = 28;
       L2 = 28;
end

if Hes == 0
    aa = L1*sqrt(abs(c2)+sqrt(abs(c4)));
    a = x + c1 - aa;
    b = x + c1 + aa;
else %Hes ==1)
    aa = sqrt(abs(c2));
    a = x + c1 - L2*aa;
    b = x + c1 + L2*aa;
end


%Grid_i = repmat((0:Ngrid-1)',1,Nstrike);    % Grid index

vk_p = @(x) calcvkp(x,b,a,strike);          % coefficients for put

w=sigma^2*T;








% fk_i = feval(rnCHF, pi*Grid_i./(b-a));
% fk_i = double(2./(b - a) .* real( fk_i ...
%     .* exp(1i .* pi .* Grid_i .* x ./ (b - a)) ...
%     .* exp(-1i .* (pi .* Grid_i .* a ./ (b - a))) ));
% 
% 
% Vk = vk_p(Grid_i);
% y = double(exp(-r*T) ...
%         .* (sum(fk_i .* Vk) - 0.5 * (fk_i(1,:) .* Vk(1,:)) ))';
% %end

if cp  == 1    % European call price using pu-call-parity
    y = y + S_0  - strike * exp(-r * T);
end

end