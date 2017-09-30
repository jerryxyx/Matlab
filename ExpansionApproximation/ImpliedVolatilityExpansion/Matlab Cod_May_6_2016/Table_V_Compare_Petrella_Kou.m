function [ ] = Table_V_Compare_Petrella_Kou( )
% This Script file calculates the prices and greeks(delta,gamma)
% of a discretely monitored barrier option for the following models(Black/Scholes,
% Variance Gamma, CGMY, Heston, Kou)
%
%
% Authors:  Guanghua Lian at Uni of South Australia
% Date: Aug 22, 2014
%
% Please send comments, suggestions, bugs, code etc. to
% ligu151@163.com; or andy.lian@unisa.edu.au
%
%
% Since this piece of code is distributed via the mathworks file-exchange
% it is covered by the BSD license
%
% This code is being provided solely for information and general
% illustrative purposes. The authors will not be responsible for the
% consequences of reliance upon using the code or for numbers produced
% from using the code.

 

Control_Limit = 5;
N_FouriersTerms = 500; 
 %% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%%  petrella 2004 numerical  Table 4 
S = 100;
K = 100;
r = 0.05;
d = 0.0;
T = 1;

UpBarrier = 101;
LowerBoundary = 0.000001; 

CallPutIndicator = -1; %% -1 means put option;  +1 means call option
InOutIndicator   = -1; %% -1 means out option; 1 means in-option
NumObservations  = 160 ; 
%% %%%%%%%%%%%%%%%%

 
% 
% %  petrella 2004 numerical  Table 9 
S = 100;
K = 100;
r = 0.05;
d = 0.0;
T = 0.2;

LowerBoundary = 0.1 ;
  
UpBarrier = 101 ;


CallPutIndicator = -1; %% -1 means put option;  +1 means call option
InOutIndicator   = -1; %% -1 means out option; 1 means in-option
NumObservations  = 50 ;
%% %%%%%%%%%%%%%%%%
 
 
 
 
 
 
 
[c p] = blsprice(S, K, r,T, 0.3);

OptionParameters = [ S, K, T, r, d, CallPutIndicator, InOutIndicator, LowerBoundary, UpBarrier, NumObservations  ];

%% %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%
 


%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 1: the upper and lower limit of the expansion;

% sigma = 0.2;
% mean = x0 + (r- div -1/2*sigma^2) *T;
% std  = sigma * sqrt(T);
% a = mean  - Control_Limit * std ;
% b = mean  + Control_Limit*std ;


%% %%%%%%%%%%%%%%%%%%%
%% CGMY
modelCGMY.ID = 'CGMY';          %model identifyer
modelCGMY.params = [.1;        %C
    .75;          %G
    1.0;          %M
    .25];         %Y

%% Parameters of the Feng (2007) Mathematical Finance paper
modelCGMY.params = [4;        %C
    50;          %G
    60;          %M
    .7];         %Y

%% %%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%
%% NIG
modelNIG.ID = 'NIG';          %model identifyer
params = [ 15, -5, 0.5];
modelNIG.params = params;
%% %%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%
%% NIG
modelNIG.ID = 'VarianceGamma';          %model identifyer
params = [ 15, -5, 0.5];
modelNIG.params = params;
%% %%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%
%% Merton
modelMerton.ID = 'Merton';          %model identifyer
%params = [ 0.1, -0.05, 0.086, 3];  %% Table 2  Feng Linetsky 2007

% params = [ 0.212, -0.01, 0.141, 2.24];     %% Table 9 of Petrella and Kou (2004)
  params = [ 0.353, -0.01, 0.141, 6.22];   %% Table 9 of Petrella and Kou (2004)
%  params = [ 0.2, 0.045, 0.3, 2];            %% Table 13 of Petrella and Kou (2004)
modelMerton.params = params;
%% %%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%
%% Black-Scholes
modelBlackScholes.ID = 'BlackScholes';          %model identifyer
params = [ 0.3 ];
modelBlackScholes.params = params;
%% %%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kou Jump model  %%
%% Note the parameter:
modelKou.ID = 'Kou';
% params = [ 0.1, 3, 0.3, 40, 12 ];  %% Table 2  Feng Linetsky 2007


rho = 0.9;
params = [ 0.0001,  48, 0.5, 5, 5 ];  %% Table 2 parameters  the parameters in Banking Finance 2013
% ratio = 0.9^2 /(1-0.9^2); lambda = params(2);  p = params(3); eta1= params(4); eta2= params(5);
%ratio = 0.1^2 /(1-0.1^2); lambda = params(2);  p = params(3); eta1= params(4); eta2= params(5);
ratio = rho^2 /(1-rho^2); lambda = params(2);  p = params(3); eta1= params(4); eta2= params(5);
params(1) =   sqrt(  ratio * 2* lambda *(   p/eta1^2 + (1-p)/eta2^2 )  );
%  
 params(1) = 0.0;

% params = [  1.011512676508229, 0, 0.5, 10, 10 ];  %% Table 2 parameters  the parameters in Banking Finance 2013
modelKou.params = params;
%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%






%% array containing all defined models

%  model = modelBlackScholes;
 model = modelMerton;
%  model = modelKou;
%model = modelCGMY;
%model = modelNIG;



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%



VanillaCallOption =     CosineMethodCallPricingFFT( 1000, 30 , model, S, K, T, r, d);  
VanillaPutOption  =     VanillaCallOption + K * exp(-r*T) - S*exp(-d*T);

bool_Comp_Greeks = false ;
tic
  [ resultBarrierOutOptions,  Delta, Gamma ] =  Function_BarrierOption_August2014(  model, OptionParameters, N_FouriersTerms, Control_Limit , bool_Comp_Greeks   ) 
%   resultBarrierInOptions = VanillaPutOption - resultBarrierOutOptions;
toc

 



% [CallDelta, PutDelta] = blsdelta( S0, K, r, T, 0.3, d);

 
format long


% ResultData = [ resultBarrierOutOptions,  Delta, Gamma ]
% save  UpOutPut.txt -ascii ResultData ;



end

