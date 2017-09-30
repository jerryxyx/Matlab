function [ Call_option ] = BarrierOption_August2014_BSM( n )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Barrier option
%% Black Scholes model here, but it can be extended to other Levy process
%% August 2 2014
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

[  S0, strike, r, div, sigma, T, l, u, n_notused ] = inputParameters( )
x0 =   0; 
% div = 0;

lowerBoundary = l ;
upBoundary = u ;

lB = log( lowerBoundary /S0  ); 
uB = log( upBoundary / S0); 
strike_X_Space = log( strike / S0 ); 


tic 

DeltaT  = T / n; 
%% ------------------------------------------------------------------------------------------------------------------- %%
%% Method 0 %%
C0 = blsprice( S0, strike, r, T, sigma)

%% ------------------------------------------------------------------------------------------------------------------- %%







%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%  Below I am using my COS method for Black-Scholes or other Levy process



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 1: the upper and lower limit of the expansion;
N_FouriersTerms = 1000 ;

mean = x0 + (r-1/2*sigma^2) *T;
std  = sigma * sqrt(T);

Control_Limit = 9;

a = mean  - Control_Limit * std ;
b = mean  + Control_Limit*std ;

N_vec = 0 : N_FouriersTerms ;
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 2: COS-expansion of the payoff function;
%% For a call option, we need to the following boundaries:
%% 1. [ log(K), infinity )
%% 2. [ a, b ]
%% 3. [ log( lowerBoundary),  log( upBoundary) ]

l = max( [ a, strike_X_Space,  lB ]);
h = min( [b,  uB ] );
%% if d < c; then we set d=c;
h = max( l, h);  

[ payoff_basis_integration ] =  f_k_expansion(   a, b, l, h, strike,  S0,     N_FouriersTerms    ) ;
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


 [ Akj_Matrix ] = A_Maxtrix_LevyProcesses(  a, b, l,  h,   DeltaT,  N_FouriersTerms  );

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







 