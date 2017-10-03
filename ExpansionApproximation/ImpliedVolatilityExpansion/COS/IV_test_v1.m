function y = IV_test_v1(S_0,strikes,T,r,q,c1,c2,c4,numGrid,isCall,isHes)
%Summary of this function goes here
% Reformat List:
%   strike -> strikes
%   cp -> isCall
%   Hes -> isHes
%   Ngrid -> numGrid
%   Nstrike -> numStrikes
% Parameters Explanation:
%   strikes: a column vector of the size Nstrike*1
%   Ngrid: the number of sample that we want to do the DFT, equals to 2^P
%   isHes: 1 for Heston model, 0 for BSM model
%   c1,c2,c4: the n-th cumulant of ln(ST/K), can be obtained using Maple.
% Example:
%>>generate_hyperparameters(50,(44:1:56)',0.1,0.01,0,0.25,0,2^6)
%>>load("hyperparamters")
%>>IV_test_v1(S_0,strikes,T,r,q,c1,c2,c4,numGrid,1,isHes)
 
x = repmat(double(log(S_0 ./ strikes))',numGrid,1);       % center
numStrikes = size(strikes,1);
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

if isHes == 0
    aa = L1*sqrt(abs(c2)+sqrt(abs(c4)));
    a = x + c1 - aa;
    b = x + c1 + aa;
else %isHes ==1)
    aa = sqrt(abs(c2));
    a = x + c1 - L2*aa;
    b = x + c1 + L2*aa;
end


Grid_i = repmat((0:numGrid-1)',1,numStrikes);    % Grid index

vk_p = @(x) calcvkp(x,b,a,strikes);          % coefficients for put

%   V_ks is a matrix indexed by grid_id(k) and strike_id(s)
V_ks = vk_p(Grid_i);







y = 0;% todo: fk_i

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

if isCall  == 1    % European call price using pu-call-parity
    y = y + S_0 * exp(-q*T)  - strikes * exp(-r*T);
end
y=V_ks;% todo: remove
end