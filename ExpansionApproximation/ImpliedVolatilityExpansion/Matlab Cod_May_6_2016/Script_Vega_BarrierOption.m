
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




[  Control_Limit, N_FouriersTerms ,  S, K, r,  d, T, LowerBoundary, UpBarrier, CallPutIndicator, InOutIndicator, NumObservations  ] = inputParameters( );
S_vector = K*(0.75:0.01:1.25);

ResultData  = zeros(  length( S_vector), 3 );

length( S_vector)
 
for i_S_vector = 1: length( S_vector)
    i_S_vector
    S = S_vector(  i_S_vector );


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
params = [ 0.1, -0.05, 0.086, 3];
params = [ 0.212, -0.01, 0.141, 2.24]; %% Table 9 of Petrella and Kou (2004)

modelMerton.params = params;
%% %%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%
%% Black-Scholes
modelBlackScholes.ID = 'BlackScholes';          %model identifyer
params = [ 0.2 ];
modelBlackScholes.params = params;
%% %%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kou Jump model  %%
%% Note the parameter:
modelKou.ID = 'Kou';
params = [ 0.1, 3, 0.3, 40, 12 ];




params = [  1.011512676508229,  48, 0.5, 5, 5 ];  %% Table 2 parameters  the parameters in Banking Finance 2013

ratio = 0.9^2 /(1-0.9^2); lambda = params(2);  p = params(3); eta1= params(4); eta2= params(5);
params(1) =   sqrt(  ratio * 2* lambda *(   p/eta1^2 + (1-p)/eta2^2 )  );
 

% params = [  1.011512676508229, 0, 0.5, 10, 10 ];  %% Table 2 parameters  the parameters in Banking Finance 2013
modelKou.params = params;
%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%






%% array containing all defined models
model = modelCGMY;
%model = modelNIG;
% model = modelMerton;
%model = modelBlackScholes;
%   model = modelKou;



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%



%result matrices
% resultBarrierOptions = zeros(length(modelVec),length(K),3);
 

VanillaCallOption = CosineMethodCallPricingFFT( 5000, 40 , model, S, K, T, r, d);  
VanillaPutOption  = VanillaCallOption + K * exp(-r*T) - S*exp(-d*T);

eps = 0.0001;

sigma = model.params(1);
% model.params(1) = sqrt(  sigma^2 + eps  );
model.params(1) =  (  sigma  + eps  );
[ resultBarrierOutOptions,  Delta, Gamma ] =  Function_BarrierOption_August2014(  model, OptionParameters, N_FouriersTerms, Control_Limit   );


model.params(1) =  (   sigma  - eps  );
[ resultBarrierOutOptions2,  Delta, Gamma ] =  Function_BarrierOption_August2014(  model, OptionParameters, N_FouriersTerms, Control_Limit   );


model.params(1) = sigma;

resultBarrierInOptions = VanillaPutOption - resultBarrierOutOptions;

% [CallDelta, PutDelta] = blsdelta( S0, K, r, T, 0.3, d);

ResultData(  i_S_vector,: ) = [ resultBarrierOutOptions,  Delta, Gamma ];

Vega(  i_S_vector ) =  (  resultBarrierOutOptions - resultBarrierOutOptions2 ) / (2*eps) ;


end

%  save  UpOutPut.txt -ascii ResultData ;
format long

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the values of the options
%  plot(  S_vector,  ResultData( :, 1)  )
%  hold on
% if strcmp(model.ID,'BlackScholes')
%     [call put] = blsprice( S_vector, K, r, T, model.params(1), d ); 
% %     plot( S_vector,  put, 'r*')
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend('m=1', 'm=2', 'm=10', 'm=52', 'm=104', 'm=252')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the Delta %%%%%%%%%%%%%%%%%%
% figure
%  plot(  S_vector,  ResultData( :, 2)  )
%  hold on
% if strcmp(model.ID,'BlackScholes')
% %     [call put] = blsdelta( S_vector, K, r, T, model.params(1), d ); 
% %     plot( S_vector,  put, 'r*')
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the Gamma %%%%%%%%%%%%%%%%%%
% figure
%  plot(  S_vector,  ResultData( :, 3)  )
 
 plot(  S_vector,  Vega   )
 
 hold on
% if strcmp(model.ID,'BlackScholes')
%    [callputblsgamma] = blsgamma( S_vector, K, r, T, model.params(1), d  ); 
%     plot( S_vector,  callputblsgamma, 'r*');
% end
% axis([75 125 -0.02 0.08]);


xlabel('Prices of the underlying stock')
ylabel('Prices of the double-barrier put option')
title('Prices of a double-barrier put option with different monitoring frequencies')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(model.ID,'BlackScholes')
    put  = blsvega( S_vector, K, r, T, model.params(1), d ); 
%     plot( S_vector,  put, 'r*')
end






