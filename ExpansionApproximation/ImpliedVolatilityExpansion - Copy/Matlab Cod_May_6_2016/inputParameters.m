

function [  Control_Limit, N_FouriersTerms ,  S, K, r,  d, T, LowerBoundary, UpBarrier, CallPutIndicator, InOutIndicator, NumObservations  ] = inputParameters( )
%% l is the barrier for the down-and-out option
%% n is the total sampling times 


%% %%%%%%%%%%%%%%%%
%% Control variables
Control_Limit =  3 ;
N_FouriersTerms = 800 ;
%% %%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fusai 2005 Finance and Stochastic;  2010 Mathematical Finance

% S = 100; 
% K = 100;
% r =0.1;
% d =0; 
% sigma =0.3;
% T = 0.2;
% 
% LowerBoundary = 99 ;
% UpBarrier = 70000 ;
% CallPutIndicator =  1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  =  25; 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %
% % EUR/AUD FX in Susuan paper Table 5
% S0 = 1.6411;
% r = 0.0838;
% q =0.0503;
% 
% 
% v0 = 0.0117;
% rho = 0.1558;
% vbar = 0.0138;
% lambda = 2.6032;
% eta = 0.3802;  %% sigma_V
% 
% 
% K = 1.477;
% LBarrier = 1.559;
% 
% T = 1.0082;
%  
% N_Observation =5;








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Feng Mathematical Finance 2007 parameter
 

S = 100; 
K = 100;
r =0.05;
d = 0.02; 
sigma =0.0;
T = 1;

LowerBoundary = 80 ;
UpBarrier = 120 ;
CallPutIndicator =  - 1; %% -1 means put option;  +1 means call option
InOutIndicator   = -1; %% -1 means out option; 1 means in-option
NumObservations  = 252 ; 

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time-dependent Blackscholes Journal of Risk
%  
% S = 100;
% K = 101;
% r = 0.05;
% d = 0.02;
% T = 0.5;
% 
% LowerBoundary = 0.1 ;
% UpBarrier = 101;
% CallPutIndicator =  -1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  = 10000 ; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Time-dependent Blackscholes Journal of Risk
 
% S = 50 ;
% K = 50 ;
% r = 0.1;
% d = 0.05;
% d = 0.00;
% T = 1;
% 
% LowerBoundary = 0.1 ;
% LowerBoundary = 40 ;
% 
% UpBarrier = 7000 ;
% CallPutIndicator =  1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  = 10 ; 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Banking & Finance %% %%%%%%%%%%%%%%%%%%%%%
% S0= 100; 
% strike = 100;
% r = 0.05;
% div = 0.0;
% sigma=0.353;
% T = 0.2;
%  
% % l= 95;
% n = 50;
% 
% % l = 105;
% l=1; 
% u = 113 ;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table 1 Journal of Banking and Finance 2013
% S0= 90; 
% strike = 96;
% r = 0.1;
% div = 0.0;
% sigma=0.3;
% T = 0.2;
%  
% % l= 95;
% n = 5;
% 
% % l = 105;
% l= 1; 
% u = 92 ;



%% Table 2 Journal of Banking and Finance 2013
% S0= 90; 
% strike = 96;
% r = 0.1;
% div = 0.0;
% sigma = 1;
% T = 0.2;
%  
% % l= 95;
% n = 5;
% 
% % l = 105;
% l= 1; 
%  l = 96 ;
%  u = 9600 ;
%  
 
 
 
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
 





%% option parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%
%  Table 2 Journal of Banking and Finance 2013
% S = 90;
% K = 96;
% r = 0.1;
% d = 0.0;
% T = 0.2;
% 
% LowerBoundary = 1 ;
% UpBarrier = 2000;
% 
% % LowerBoundary = 96 ;
% UpBarrier = 96;
% 
% CallPutIndicator = -1; %% -1 means put option;  +1 means call option
% InOutIndicator   = 1; %% -1 means out option; 1 means in-option
% NumObservations  = 5 ;
%% End %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%%  petrella 2004 numerical  Table 4 
% S = 100;
% K = 100;
% r = 0.05;
% d = 0.0;
% T = 1;
% 
% LowerBoundary = 0.1 ;
% UpBarrier = 2000;
% 
% % LowerBoundary = 96 ;
% UpBarrier = 101;
% UpBarrier = 101;
% 
% CallPutIndicator = -1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  = 5 ; 
%% %%%%%%%%%%%%%%%%






%%  petrella 2004 numerical  Table 7 
% S = 100;
% K = 100;
% r = 0.05;
% d = 0.0;
% T = 0.5;
% 
% LowerBoundary = 0.1 ;
% UpBarrier = 2000;
%  
% UpBarrier = 101;
% UpBarrier = 100.05;
% 
% CallPutIndicator = -1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  = 50 ; 
%% %%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %  petrella 2004 numerical  Table 9 
% S = 100;
% K = 100;
% r = 0.05;
% d = 0.0;
% T = 0.2;
% 
% LowerBoundary = 0.1 ;
% UpBarrier = 2000;
%  
% UpBarrier = 113 ;
% 
% 
% CallPutIndicator = -1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  = 50 ;
%% %%%%%%%%%%%%%%%%

%%  petrella 2004 numerical  Table 13
% S = 100;
% K = 100;
% r = 0.1;
% d = 0.0;
% T = 0.2;
% 
% LowerBoundary = 99 ;
% UpBarrier = 2000;
%  
% % UpBarrier = 105;
% 
% 
% CallPutIndicator =  1; %% -1 means put option;  +1 means call option
% InOutIndicator   = -1; %% -1 means out option; 1 means in-option
% NumObservations  = 5 ; 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















 
