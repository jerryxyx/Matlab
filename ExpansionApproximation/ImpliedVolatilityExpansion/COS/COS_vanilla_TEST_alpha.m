function y = COS_vanilla_TEST_alpha(rnCHF,r,T,S_0,cp,strike,Ngrid,Hes,alpha)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

Nstrike = size(strike,1);

x = repmat(double(log(S_0 ./ strike))',Ngrid,1);       % center


a = x - alpha;
b = x + alpha;

Grid_i = repmat((0:Ngrid-1)',1,Nstrike);    % Grid index

vk_p = @(x) calcvkp(x,b,a,strike);          % coefficients for put

fk_i = feval(rnCHF, pi*Grid_i./(b-a));
fk_i = double(2./(b - a) .* real( fk_i ...
    .* exp(1i .* pi .* Grid_i .* x ./ (b - a)) ...
    .* exp(-1i .* (pi .* Grid_i .* a ./ (b - a))) ));


Vk = vk_p(Grid_i);
y = double(exp(-r*T) ...
        .* (sum(fk_i .* Vk) - 0.5 * (fk_i(1,:) .* Vk(1,:)) ))';
%end

if cp  == 1    % European call price using pu-call-parity
    y = y + S_0  - strike * exp(-r * T);
end

end