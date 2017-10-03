function [ Call_option ] = BarrierOption_August2014_CGMY( n )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Barrier option
%% Black Scholes model here, but it can be extended to other Levy process
%% August 2 2014
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
[  S0, strike, r, div, sigma, T, l, u, n_NOTUSED ] = inputParameters( );

% n = 252;

x0   =   0; 
% div = 0;

lowerBoundary = l ;
upBoundary = u;

% lowerBoundary = 0.1 ;
% upBoundary = 1000;


lB = log( lowerBoundary /S0  ); 
uB = log( upBoundary / S0); 
strike_X_Space = log( strike / S0 ); 


tic 

DeltaT  = T / n; 
%% ------------------------------------------------------------------------------------------------------------------- %%
%% Method 0 %%
C0 = blsprice( S0, strike, r, T, sigma, div)

%% ------------------------------------------------------------------------------------------------------------------- %%

 

%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%  Below I am using my COS method for Black-Scholes or other Levy process



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 1: the upper and lower limit of the expansion;
N_FouriersTerms = 5000 ;
sigma = 0.2;
mean = x0 + (r- div -1/2*sigma^2) *T;
std  = sigma * sqrt(T);



Control_Limit = 50;

a = mean  - Control_Limit * std ;
b = mean  + Control_Limit*std ;

sigma = 0;
 

%% %%%%%%%%%%%%%%%%%%% 
%% CGMY
modelCGMY.ID = 'CGMY';          %model identifyer
modelCGMY.params = [.1;        %C 
                      .75;          %G
                      1.0;          %M
                      .25];         %Y
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
modelMerton.ID = 'BlackScholes';          %model identifyer
  params = [ 0.1, -0.05, 0.086, 3];
modelMerton.params = params;
%% %%%%%%%%%%%%%%%%%%%   


%% %%%%%%%%%%%%%%%%%%% 
%% Black-Scholes
modelBlackScholes.ID = 'BlackScholes';          %model identifyer
  params = [ 0.5 ];
modelBlackScholes.params = params;
%% %%%%%%%%%%%%%%%%%%%   


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kou Jump model  %%     
%% Note the parameter: 
modelKou.ID = 'Kou';  
        params = [ 0.1, 3, 0.3, 40, 12 ];
modelKou.params = params; 
%% %%%%%%%%%%%%%%%%%%%   
%% %%%%%%%%%%%%%%%%%%%   

        
%array containing all defined models
model = modelCGMY;
model = modelNIG;
model = modelMerton;
model = modelBlackScholes;
model = modelBlackScholes;



[c1, c2, c4] = getCumulants(T,r,div,model);

    % truncation range
    a = c1 - Control_Limit*sqrt(abs(c2) + sqrt(abs(c4)));
    b = 2*c1 - a;

%  phi = exp(feval(@CharacteristicFunctionLib, model,tmp,T,r,d,model.params));
 
N_vec = 0 : N_FouriersTerms ;
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 2: COS-expansion of the payoff function;
%% For a call option, we need to the following boundaries:
%% 1. [ log(K), infinity )
%% 2. [ a, b ]
%% 3. [ log( lowerBoundary),  log( upBoundary) ]

% l = max( [ a, strike_X_Space,  lB ]);
% h = min( [b,  uB ] );
% % if d < c; then we set d=c;
% h = max( l, h);  
% [ payoff_basis_integration ] =  f_k_expansion(   a, b, l, h, strike,  S0,     N_FouriersTerms    ) ;

%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% For a put option, we need to the following boundaries:
%% 1. [ -infinity, log(K)  )
%% 2. [ a, b ]
%% 3. [ log( lowerBoundary),  log( upBoundary) ]


l = max( [ a,   lB  ]);
h = min( [b, strike_X_Space,   uB  ] );
%% if d < c; then we set d=c;
h = max( l, h);  
[ payoff_basis_integration ] =  - f_k_expansion(   a, b, l, h, strike,  S0,     N_FouriersTerms    ) ;

%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%





%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 3: COS-expansion of the probability density function;
%% We need to the following boundaries:
%% 2. [ a, b ]
%% 3. [ log( lowerBoundary),  log( upBoundary) ]

l = max( [ a,   lB ]);
h = min( [b,  uB ] );
%% if d < c; then we set d=c;
h = max( l, h);   

% mu_deltaT    = ( r- div - 0.5*sigma^2 ) * T /n ;
% sigma_deltaT = (sigma) * sqrt(T/n);

  [  A_j ] = A_j_Vector(   a, b, x0,  DeltaT,   N_FouriersTerms    );
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 4: COS-expansion of the probability density function;
%% So we can compute a vector to get a vector for the final step computation

% l = max( [ a,  log( lowerBoundary) ] )  ;
% u = min( [ b, log( upBoundary)] )  ;
% %% if d < c; then we set d=c;
% u = max( l, u); 
% 
% mu_deltaT    = ( r- div - 0.5*sigma^2 ) * T /n ;
% sigma_deltaT = (sigma) * sqrt(T/n);
 
tic
 [ Akj_Matrix ] = A_Maxtrix_LevyProcesses(  a, b, l,  h,   DeltaT,  N_FouriersTerms  );
toc
%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 5: Find out the option price by:
%% C = exp( - r*T) * [  coefficient_vec * COS_jk_Integrate^(n-1) * payoff_basis_integration']
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%

C2  = exp( -r*T ) * A_j'  * Akj_Matrix^(n-1) * payoff_basis_integration  


%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%


%% End of my COS method for Black-Scholes or other Levy process
%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%

toc

end







 