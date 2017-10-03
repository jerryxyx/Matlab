function [ OptionValue, Delta, Gamma ] = Function_BarrierOption_August2014(  model, OptionParameters, N_FouriersTerms, Control_Limit, bool_Comp_Greeks    )

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Barrier option
%% Black Scholes model here, but it can be extended to other Levy process
%% August 2 2014
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% OptionParameters = [ S, K, T, r, d, CallPutIndicator, InOutIndicator, LowerBoundary, UpBarrier, NumObservations  ];

S0 = OptionParameters(1);
K = OptionParameters(2);
T = OptionParameters(3);
r = OptionParameters(4);
d = OptionParameters(5);

CallPutIndicator    = OptionParameters(6);
InOutIndicator      = OptionParameters(7);
LowerBoundary       = OptionParameters(8);
UpBarrier           = OptionParameters(9);
NumObservations     = OptionParameters(10);


[c1, c2, c4] = getCumulants(T,r,d,model);
% truncation range
a = c1 + log(S0) - Control_Limit*sqrt(abs(c2) + sqrt(abs(c4)));
b = c1 + log(S0) + Control_Limit*sqrt(abs(c2) + sqrt(abs(c4)));
%  phi = exp(feval(@CharacteristicFunctionLib, model,tmp,T,r,d,model.params));



x0   =   log(S0) ;
%  x0   =  0;
S0 = 1;

lB = log( LowerBoundary /S0  );
uB = log( UpBarrier / S0);
strike_X_Space = log( K / S0 );
DeltaT  = T / NumObservations;


% N_vec = 0 : N_FouriersTerms ;

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




%% If it is a put option
if CallPutIndicator < 0
    l = max( [ a,   lB  ]);
    h = min( [b, strike_X_Space,   uB  ] );
    %% if d < c; then we set d=c;
    h = max( l, h );
    [ payoff_basis_integration ] =   f_k_expansion(   a, b, l, h, K,  S0,     N_FouriersTerms    ) ;
    payoff_basis_integration = - payoff_basis_integration;
    
elseif CallPutIndicator > 0  %% if it is a call option
    
    l = max( [ a,  strike_X_Space,  lB  ]);
    h = min( [b, uB  ] );
    %% if d < c; then we set d=c;
    h = max( l, h);
    
    [ payoff_basis_integration ] =   f_k_expansion(   a, b, l, h, K,  S0,     N_FouriersTerms    ) ;
end




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

[  A_j, Delta_A_j] = A_j_Vector(   a, b, x0,  DeltaT,   N_FouriersTerms,   model, OptionParameters, bool_Comp_Greeks   );
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bool_Comp_Greeks == true
    Gamma_A_j =   - A_j   .* ( (0:N_FouriersTerms).^2 )' *  (pi/(b-a))^2;
end

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




[ Akj_Matrix ] = A_Maxtrix_LevyProcesses(  a, b, l,  h,   DeltaT,  N_FouriersTerms,  model, OptionParameters  );

% tic
% [ Akj_Matrix_SecondPart ] = A_Maxtrix_SecondPart_LevyProcesses(  a, b, l,  h,  N_FouriersTerms   );
% [ Akj_Matrix ] = A_Maxtrix_LevyProcesses_version2(  a, b, l,  h,   DeltaT,  N_FouriersTerms,   model, OptionParameters, Akj_Matrix_SecondPart   );
% toc


%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
%% Part 5: Find out the option price by:
%% C = exp( - r*T) * [  coefficient_vec * COS_jk_Integrate^(n-1) * payoff_basis_integration']
%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%


% FirstPowerPart = Akj_Matrix^( NumObservations-1) * payoff_basis_integration;


if NumObservations< 4000
    
    FirstPowerPart = payoff_basis_integration;
    
    for index =1 : NumObservations-1
        FirstPowerPart = Akj_Matrix* FirstPowerPart;
    end
    
    
elseif  NumObservations >= 4000
    % tic; FirstPowerPart = Akj_Matrix^( NumObservations-1) * payoff_basis_integration;  toc
    
    FirstPowerPart= mpower2(Akj_Matrix,  NumObservations-1 ) * payoff_basis_integration;
    
end


%
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% A new method to compute the Akj_Matrix^( NumObservations-1)
%  BinNum = dec2bin(  NumObservations-1 ) .' - '0' ;
%
%  PositionPositive = length( BinNum) - find( BinNum >0 ) ;
%
% sizeOfMatrix = size( Akj_Matrix );
%
%  A = ones( sizeOfMatrix(1), sizeOfMatrix(2) , length( find( BinNum >0 ) ));
%
%  ithPosition = 1;
%
%  if find(0==PositionPositive)
%      A(:,:,1) = Akj_Matrix ;
%      ithPosition = ithPosition +1 ;
%  end
%
%
%  NumOfPosition =  length( find( BinNum >0 ) );
%
% tic
%  A(:,:, NumOfPosition )  = Akj_Matrix;
%
%  for i = 1: PositionPositive(1)
%     A(:,:, NumOfPosition )  = A(:,:, NumOfPosition ) ^2;
%
%     if  find( PositionPositive == i )
%         A(:,:,ithPosition ) = A(:,:, NumOfPosition )  ;
%         ithPosition = ithPosition +1 ;
%     end
%
%  end
%  toc
%
%  tic
%
%  FirstPowerPart = payoff_basis_integration;
% for index =1 : NumOfPosition
%     FirstPowerPart =  A(:,:, index )* FirstPowerPart;
% end
% toc
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



OptionValue  = exp( -r*T ) * A_j'  *  FirstPowerPart;

Delta =0;
Gamma = 0;
if bool_Comp_Greeks == true
    Delta = exp( -r*T ) * Delta_A_j'  * FirstPowerPart ;
    Gamma = exp( -r*T ) / exp(2*x0) * ( Gamma_A_j' - Delta_A_j' * exp(x0) ) * FirstPowerPart;
end

%% %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%


%% End of my COS method for Black-Scholes or other Levy process
%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%


end







