addpath('E:\Matlab Code\Pricing by Integration\RN_CHFs')

Hes = 0;
r   = .05;
T   = 1/12;
S_0 = 50;
F_0 = S_0*exp(r*T);

%STRIKES
cp = 1;  %cp ==1 is price calls
W = (44:1:56);
strike = W';
K = size(strike,1);

%GRID
P = 6;
Ngrid = 2 ^ P;

if T>=.1
    L1 = 10;
    L2 = 12; %for Hes
else
    L1 = 15;
    L2 = 18;
end


%-----------------RN CHF/PARAMETERS-----------------------------------%
% 
% %%%-----Merton CHF: 
% sigma = .3;
% muj = .1;
% sigmaj = .2;
% lam = .5;
% rnCHF = @(u)cf_RN_MJD( u, r, T, sigma, muj, sigmaj , lam);

%RNmu = r - .5*sigma^2-lam*(exp(muj - .5*sigmaj^2)-1);
%c1 = T*(RNmu+lam*muj);
%c2 = T*lam*(sigma^2/lam + muj^2 +sigmaj^2);
%c4 = T*lam*(muj^4 + 6*sigmaj^2*muj^2+3*sigmaj^4);

%%%------------------------------------------------------------------------

% %%%%-----NIG CHF
% alpha = 3;
% beta = 0;
% delta = 5;
% rnCHF = @(u)cf_RN_NIG( u,r,T,alpha,beta,delta);
% 
% asq = alpha^2; bsq = beta^2;
% temp = asq - bsq;
% RNmu = r + delta*(sqrt(asq - (beta+1)^2)-sqrt(temp));
% c1 = T*(RNmu + delta*beta/sqrt(temp));
% c2 = delta*T*asq*temp^(-1.5);
% c4 = 3*delta*T*asq*(asq + 4*bsq)*(temp)^(-3.5);

%  
%%%------------------------------------------------------------------------

% %-------VG CHF  (USE T=1)
% C = 4.3;   %Default: (C,G,M) = (4.3,5.2,5.6), T=1, r=.05
% G = 52;
% MM = 50;
% rnCHF = @(u) cf_RN_VG(u,T,r,C,G,MM);

%RNmu = r - C*log(G*MM/(G*MM + (MM-G) -1));
%c1 = ;
%c2 = ;
%c4 = ;

%%%------------------------------------------------------------------------

% % %-------VG2 CHF  (USE T=1)
% sigma = .12;
% theta = -.14;
% nu = .2;
% rnCHF = @(u)cf_RN_VG2(u,T,r,sigma,nu,theta);
% 
% sig2 = .5*sigma^2;
% RNmu = r + 1/nu*log(1-nu*(theta+sig2));
% c1 = (RNmu + theta)*T;
% c2 = (sig2+nu*theta^2)*T;
% c4 = 3*(sig2^2*nu +2*theta^4*nu^3+4*sig2*theta^2*nu^2)*T;


% %%%------------------------------------------------------------------------
% %%%-------CGMY CHF  (USE T=1)
% C = 1;   %Default: (C,G,M,Y) = (4.3,5.2,5.6,.1), T=1, r=.05
% G = 5;
% MM = 5;
% Y = .7;  %(Y<=.5 for stability)
% rnCHF = @(u)cf_RN_CGMY(u,T,r,C,G,MM,Y);
% 
% RNmu = r - C*gamma(-Y)*((MM-1)^Y - MM^Y + (G+1)^Y - G^Y); 
% c1 = RNmu*T + C*T*gamma(1-Y)*(MM^(Y-1)-G^(Y-1));
% c2 = C*T *gamma(2-Y)*(MM^(Y-2)+G^(Y-2));
% c4 = C*T *gamma(4-Y)*(MM^(Y-4)+G^(Y-4));

%%%------------------------------------------------------------------------

% %%---- BSM CHF
sigmaBSM = .3;
rnCHF = @(u) cf_RN_BSM( u, r, T, sigmaBSM );

RNmu = r-.5*sigmaBSM^2;
c1 = RNmu*T;
c2 = sigmaBSM^2*T;
c4 = 0;


%%%------------------------------------------------------------------------

% %%----Heston CHF:   Base: (u_0,ubar,lam,eta,rho)=(.02,.02,.01,.2,0)
% % u_0 = .0175;   
% % ubar = .0398;
% % lam = 1.5768;
% % eta = .5751; % <=.25
% % rho = -0.5711;
% 
% u_0 = .02;   
% ubar = .01;
% lam = .01;
% eta = .2; % <=.25
% rho = -.5;
% 
% rnCHF = @(u)cf_RN_Heston(u,T,r,u_0,ubar,lam,eta,rho);
% 
% RNmu = r;
% c1 = (RNmu)*T + (1-exp(-lam*T))*(ubar - u_0)/(2*lam);
% c2 = 1/(8*lam^3)*(eta*T*lam*exp(-lam*T)*(u_0-ubar)*(8*lam*rho-4*eta)...
%     + lam*rho*eta*(1-exp(-lam*T))*(16*ubar - 8*u_0)...
%     +2*ubar*lam*T*(-4*lam*rho*eta + eta^2 + 4*lam^2)...
%     + eta^2*((ubar - 2*u_0)*exp(-2*lam*T) + ubar*(6*exp(-lam*T)-7)+2*u_0)...
%     + 8*lam^2*(u_0-ubar)*(1-exp(-lam*T)));
% c4 = 0;
% Hes = 1;

%%%------------------------------------------------------------------------

% %-----KOU CHF  (NOOOO GOOOD)
% sigma = .03;
% lam = .01;
% p_up =.9;
% eta1 = 10;
% eta2 = 10;
% rnCHF = @(u)cf_RN_KOU(u,T,r,sigma,lam,p_up,eta1,eta2);

%RNmu = ;
%c1 = ;
%c2 = ;
%c4 = ;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = repmat(double(log(S_0 ./ strike))',Ngrid,1);       % center
%yy = zeros(size(x));

if Hes == 0
    aa = L1*sqrt(abs(c2)+sqrt(abs(c4)));
    a = x + abs(c1) - aa;
    b = x + abs(c1) + aa;
    
else %Hes ==1)
    aa = sqrt(abs(c2));
    a = x + abs(c1) - L2*aa;
    b = x + abs(c1) + L2*aa;
end
tic



Grid_i = repmat((0:Ngrid-1)',1,K);    % Grid index

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


toc

% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ---- PRICE COMPARISONS
% BSM = 0;
% CM = 0;
% H = 1;   %index of where we start including values in the error comparison
% 
% RefValues = zeros(K,1);
% if BSM == 1
%     rst = (r+ sigmaBSM^2/2)*T;
%     st = sigmaBSM*sqrt(T);
%     disc = exp(-r*T);
%     for k=1:K
%         d1 = 1/st*(log(S_0/strike(k)) + rst);
%         d2 = d1 - st;
%         RefValues(k) = S_0*normcdf(d1) - normcdf(d2)*strike(k)*disc;  
%     end
% elseif CM ==1
%     CMFT_N = 2^23;  %Need 2^25 to get e-11, 2^18 to get e-07, 2^20 to get e-08
%     for k=1:K
%         RefValues(k) = call_price_CM( rnCHF, 'cubic',CMFT_N, S_0, strike(k), T, r );
%     end
% else 
%     RefValues = COS_vanilla_TEST(rnCHF,r,T,S_0,cp,W',2^19,Hes,c1,c2,c4);
% end
% %mindev = min(abs(y(H:end)-RefValues(H:end)))
% maxdev = max(abs(y(H:end)-RefValues(H:end)));
% avgdev = mean(abs(y(H:end)-RefValues(H:end)))
