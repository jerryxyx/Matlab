function []=generate_hyperparameters(S_0,strikes,T,r,q,sigmaBSM,isHes,numGrid)

% script for generating hyperparameters

%--------------------------------------------------
% one can use parameters below
%{
S0 = 100;
r = 0.1;
q = 0;
T = 0.1;
sigmaBSM = 0.25;
Hes = 0;
strikes = (44:1:56)';
Ngrid = 2^6
%}
%--------------------------------------------------
% numGrid*numStrikes
x = repmat(double(log(S_0 ./ strikes))',numGrid,1);       % center

RNmu = r-q-.5*sigmaBSM^2;

% cumulant of ln(ST/K)
% for GBM
c1 = RNmu*T;
c2 = sigmaBSM^2*T;
c4 = 0;

% truncated interval
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

% save in parameters in "hyperparameters.mat"
save("hyperparameters");
end
