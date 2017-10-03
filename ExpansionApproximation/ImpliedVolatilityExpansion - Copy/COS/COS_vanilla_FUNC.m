function y = COS_vanilla_FUNC(rnCHF,r,T,S_0,cp,strike,Ngrid,Hes,c1,c2,c4)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

Nstrike = size(strike,1);

if T>=.1
    L1 = 10;
    L2 = 12;
else
    L1 = 15;
    L2 = 18;
end

if Hes == 0
    aa = L1*sqrt(abs(c2)+sqrt(abs(c4)));
    ab = c1 - aa;  %abar
else %Hes ==1)
    aa = sqrt(abs(c2));
    ab = c1 - L2*aa;
    aa = aa*L2;
end

ax = log(S_0 ./ strike) + ab;       % column vec

gam = pi/(2*aa);
Grid1 = (1:Ngrid-1);    % Grid index
Grid0 = [0 Grid1];
mu = exp(-1i*gam*ab*Grid0);
mu = real(feval(rnCHF,gam*Grid0).*mu)';
%M = repmat(mu,Nstrike,1);
U = zeros(Nstrike,Ngrid);
E = exp(ax);

for k=1:Nstrike
   U(k,2:end) = 1/aa*((1./(1+(gam*Grid1).^2)).*(-cos(gam*ax(k)*Grid1) + E(k) + Grid1.*gam.*sin(gam*ax(k)*Grid1)) -1./(gam*Grid1).*sin(gam*ax(k)*Grid1));
   U(k,1) = 1/(2*aa)*(E(k)-1-ax(k));
end

%M = real(U.*M);
%y = exp(-r*T)*sum(M,2).*strike;
y = exp(-r*T)*U*mu.*strike;

if cp  == 1    % European call price using pu-call-parity
    y = y + S_0  - strike * exp(-r * T);
end



end


% x = repmat(double(log(S_0 ./ strike))',Ngrid,1);       % center
% 
% if T>=.1
%     L1 = 10;
%     L2 = 12;
% else
%     L1 = 15;
%     L2 = 18;
% end
% 
% %yy = zeros(size(x));
% 
% if Hes == 0
%     aa = L1*sqrt(abs(c2)+sqrt(abs(c4)));
%     a = x + c1 - aa;
%     b = x + c1 + aa;
%     vv = 2*aa;
% else %Hes ==1)
%     aa = sqrt(abs(c2));
%     a = x + c1 - L2*aa;
%     b = x + c1 + L2*aa;
%     vv = 2*aa*L2;
% end
% 
% 
% Grid_i = repmat((0:Ngrid-1)',1,Nstrike);    % Grid index
% 
% vk_p = @(x) calcvkp(x,b,a,strike);          % coefficients for put
% 
% gam = pi/vv;
% %fk_i = feval(rnCHF, pi*Grid_i./(b-a));
% fk_i = repmat(feval(rnCHF,gam*Grid_i(:,1)),1,Nstrike);
% 
% fk_i = double(2./(b - a) .* real( fk_i ...
%     .* exp(1i .* Grid_i .* x *gam) ...
%     .* exp(-1i .* (Grid_i .* a *gam)) ));
% 
% 
% Vk = vk_p(Grid_i);
% y = double(exp(-r*T) ...
%         .* (sum(fk_i .* Vk) - 0.5 * (fk_i(1,:) .* Vk(1,:)) ))';
% %end
% 
% if cp  == 1    % European call price using pu-call-parity
%     y = y + S_0  - strike * exp(-r * T);
% end

