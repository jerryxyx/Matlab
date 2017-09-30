

function [ Akj_Matrix ] = A_Maxtrix_LevyProcesses(  a, b, l,  h,   DeltaT,  N_FouriersTerms,   model, OptionParameters   )

%% The code generates the A_jkn Matrix


%% xa, xb are the boundaries of the Fourier expansion to truncate the integration;
%% c and d are the boundaries computing the integration A_j,k,n


%%
 


%% Declair the total number of terms of the Fourier expansions
%%
N = N_FouriersTerms;


vec_j = [0:N];
vec_k = [0:N];


Akj_Matrix= zeros(N,N);




psi_fun_vector = zeros(1, length(vec_j));

exp_psi_fun_vector = zeros(1, length(vec_j));

for j_index=1:length(vec_j)
    j = vec_j(j_index);
    psi_fun_vector(j_index) = psi_fun(j*pi/(b-a),   model, OptionParameters  );
end
     
exp_psi_fun_vector = exp( DeltaT * psi_fun_vector ); 

%% ksi( ) is a function in the characteristic function exponential. 
%% If k == j then ...
for k_index=1:length(vec_k)
    k = vec_k(k_index);
    Akj_Matrix(k_index, k_index) =  2/(b-a) * real(  exp( DeltaT * psi_fun_vector(k_index) ) * ( (((i) * a + (-1*i) * b) * sin(pi * k * (-l + a) / (-b + a)) ^ 2 + cos(pi * k * (-l + a) / (-b + a)) * (-b + a) * sin(pi * k * (-l + a) / (-b + a)) + ((-1*i) * a + (i) * b) * sin(pi * k * (-h + a) / (-b + a)) ^ 2 - cos(pi * k * (-h + a) / (-b + a)) * (-b + a) * sin(pi * k * (-h + a) / (-b + a)) - pi * k * (l - h)) / pi / k / 2 )   );
end
 
%% Else if k != j ...
% for k_index=1:length(vec_k)
%     k = vec_k(k_index);
%      for j_index=1:length(vec_j)
%         j = vec_j(j_index);   
%         if( j~=k)
%              Akj_Matrix(k_index, j_index)  = 2/(b-a) * real(  exp( DeltaT * psi_fun(j*pi/(b-a) ) ) * (-((i) * j * exp((i) * pi * j * (-l + a) / (-b + a)) * cos(pi * k * (-l + a) / (-b + a)) + k * exp((i) * pi * j * (-l + a) / (-b + a)) * sin(pi * k * (-l + a) / (-b + a)) + (-1*i) * j * exp((i) * pi * j * (-h + a) / (-b + a)) * cos(pi * k * (-h + a) / (-b + a)) - k * exp((i) * pi * j * (-h + a) / (-b + a)) * sin(pi * k * (-h + a) / (-b + a))) * (-b + a) / pi / (j ^ 2 - k ^ 2) ) );
%         end
%         
%     end
% end

%% New version, updated on Aug 19


COS_k_vector_l = zeros(1, length(vec_k));
SIN_k_vector_l = zeros(1, length(vec_k));
COS_k_vector_l = cos(pi * vec_k * (- l + a) / (-b + a));
SIN_k_vector_l  = sin(pi * vec_k * (- l + a) / (-b + a));

COS_k_vector_u = zeros(1, length(vec_k));
SIN_k_vector_u = zeros(1, length(vec_k));
COS_k_vector_u = cos(pi * vec_k * (- h + a) / (-b + a));
SIN_k_vector_u = sin(pi * vec_k * (- h + a) / (-b + a));

EXP_j_vector_u = zeros(1, length(vec_j));
EXP_j_vector_u = exp( 1i * pi * vec_j * (- h + a) / (-b + a));
 
EXP_j_vector_l = zeros(1, length(vec_j));
EXP_j_vector_l = exp( 1i * pi * vec_j * (- l + a) / (-b + a));


for k_index=1:length(vec_k)
    k = vec_k(k_index);
     for j_index=1:length(vec_j)
        j = vec_j(j_index);   
        if( j~=k)
%             Akj_Matrix(k_index, j_index)  = 2/(b-a) * real(  exp( DeltaT * psi_fun_vector(j_index)  ) * (-((i) * j * exp((i) * pi * j * (-l + a) / (-b + a)) * cos(pi * k * (-l + a) / (-b + a)) + k * exp((i) * pi * j * (-l + a) / (-b + a)) * sin(pi * k * (-l + a) / (-b + a)) + (-1*i) * j * exp((i) * pi * j * (-h + a) / (-b + a)) * cos(pi * k * (-h + a) / (-b + a)) - k * exp((i) * pi * j * (-h + a) / (-b + a)) * sin(pi * k * (-h + a) / (-b + a))) * (-b + a) / pi / (j ^ 2 - k ^ 2) ) );
            Akj_Matrix(k_index, j_index)  = 2 / pi / (j ^ 2 - k ^ 2)  * real( exp_psi_fun_vector(j_index) * (  (1i*j* COS_k_vector_l(k_index) + k*SIN_k_vector_l(k_index) )* EXP_j_vector_l(j_index)  - ( 1i*j* COS_k_vector_u(k_index) + k*SIN_k_vector_u(k_index)  )* EXP_j_vector_u(j_index) )   ) ;  
        end
        
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%  Akj_Matrix(1, 1) =  2/(b-a) * real(  exp( DeltaT * psi_fun( k*pi/(b-a) ) ) * (h-l ) );
  Akj_Matrix(1, 1) =    2/(b-a) * ( h-l ) ;
 %% Because the A_j is the coefficient of the Fourier series, A_0 should be divided by 2; 
 Akj_Matrix( :, 1) =  Akj_Matrix( :, 1)  / 2 ;  
 
%  
% [ S0, strike, r, sigma, T, lllll, n ] = inputParameters( );
% A2kj_Matrix= zeros(N,N);
% A3kj_Matrix= zeros(N,N);
% 
% %     densityNormal= @(x, z ) (  1./ ( sqrt( 2.*pi* sigma.^2 .* DeltaT ) ) .*  (  exp(  - ( x- ( z + ( r - sigma.^2/2 ) * DeltaT )  ).^2 ./ (2 * sigma.^2 .* DeltaT )   )  )    );
%     
%     mu_NEW = ( r - sigma.^2/2 )* DeltaT;
%     mu_deltaT = ( r - sigma.^2/2 )* DeltaT;
%     sigma_deltaT =  sigma * sqrt(DeltaT) ;
    
%  
%  for k_index=1:length(vec_k)
%     k = vec_k(k_index);
%      for j_index=1:length(vec_j)
%         j = vec_j(j_index);   
% %         f= @(x, z) ( 2 ./ (b-a)  .* ( densityNormal( x , z ) ) .* cos( j .* pi .* ( x - a) ./ (b - a)) .*  cos( k .* pi .* ( z - a) ./ (b - a)) );
% %         A2kj_Matrix(k_index, j_index) =  integral2(  f, a, b, l, h ) ;  
%         if( k~= j )
% %         A3kj_Matrix(k_index, j_index) =  2 ./ (b-a)   * exp(-j ^ 2 * pi ^ 2 * sigma ^ 2 * DeltaT / (-b + a) ^ 2 / 0.2e1) * (-b + a) / pi / (j - k) / (j + k) * (-sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j + 0.2e1 * k * a + 0.2e1 * l * j - 0.2e1 * l * k) / (-b + a) / 0.2e1) * j - sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j + 0.2e1 * k * a + 0.2e1 * l * j - 0.2e1 * l * k) / (-b + a) / 0.2e1) * k - sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j - 0.2e1 * k * a + 0.2e1 * l * j + 0.2e1 * l * k) / (-b + a) / 0.2e1) * j + sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j - 0.2e1 * k * a + 0.2e1 * l * j + 0.2e1 * l * k) / (-b + a) / 0.2e1) * k + sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j + 0.2e1 * k * a + 0.2e1 * h * j - 0.2e1 * h * k) / (-b + a) / 0.2e1) * j + sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j + 0.2e1 * k * a + 0.2e1 * h * j - 0.2e1 * h * k) / (-b + a) / 0.2e1) * k + sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j - 0.2e1 * k * a + 0.2e1 * h * j + 0.2e1 * h * k) / (-b + a) / 0.2e1) * j - sin(pi * (-DeltaT * j * sigma ^ 2 + 0.2e1 * DeltaT * j * r - 0.2e1 * a * j - 0.2e1 * k * a + 0.2e1 * h * j + 0.2e1 * h * k) / (-b + a) / 0.2e1) * k) / 0.2e1;
%         
% %         A3kj_Matrix(k_index, j_index) =  2  * exp(-j ^ 2 * pi ^ 2 * sigma ^ 2 * DeltaT / (-b + a) ^ 2 / 0.2e1) / pi / (j - k) / (j + k) * ...
% %                                             (  (j+k) * sin( pi*( (l-h)*(k-j) ) /2/(b-a) ) * cos( (  pi* ( (h+l -2*a)*(k-j) - 2*j* mu_NEW ) )  /(2*(b-a) )) + ...
% %                                                 (j-k) * sin( pi*( (h-l)*(k+j) ) /2/(b-a) ) * cos( (  pi* ( (h+l -2*a)*(k+j) + 2*j* mu_NEW ) )  /(2*(b-a) ) ))  ;
% 
%         c= l;
%         d= h;
%         A3kj_Matrix(k_index, j_index) = (2/(b-a)) *  exp( -j.^ 2 * pi ^ 2 * sigma_deltaT ^ 2 / (-b + a) ^ 2 / 0.2e1 )...
%             *(-b + a) * (j * sin(pi * (-c * j + k * c + j * a - j * mu_deltaT - k * a) / (-b + a)) + sin(pi * (-c * j + k * c + j * a - j * mu_deltaT - k * a) / (-b + a)) * k + j * sin(pi * (-c * j - k * c + j * a - j * mu_deltaT + k * a) / (-b + a)) - k * sin(pi * (-c * j - k * c + j * a - j * mu_deltaT + k * a) / (-b + a)) - sin(pi * (-d * j + d * k + j * a - j * mu_deltaT - k * a) / (-b + a)) * j - sin(pi * (-d * j + d * k + j * a - j * mu_deltaT - k * a) / (-b + a)) * k - sin(pi * (-d * j - d * k + j * a - j * mu_deltaT + k * a) / (-b + a)) * j + sin(pi * (-d * j - d * k + j * a - j * mu_deltaT + k * a) / (-b + a)) * k) / pi / (j - k) / (j + k) / 0.2e1;    
%         end
%      end 
%  end

 
%   A2kj_Matrix( :, 1) =  A2kj_Matrix(  :, 1)  / 2 ;  



    
%   for k_index=1:length(vec_k)
%     j = vec_k(k_index);
%     A3kj_Matrix(k_index, k_index) =  2/(b-a) *  ( exp(-j ^ 2 * pi ^ 2 * sigma ^ 2 * DeltaT / (-b + a) ^ 2 / 0.2e1) * (0.2e1 * cos(pi * j * (-sigma ^ 2 + (2 * r)) * DeltaT / (-b + a) / 0.2e1) * h * pi * j - 0.2e1 * cos(pi * j * (-sigma ^ 2 + (2 * r)) * DeltaT / (-b + a) / 0.2e1) * l * pi * j - sin(pi * j * (-DeltaT * sigma ^ 2 + 0.2e1 * DeltaT * r - 0.4e1 * a + 0.4e1 * l) / (-b + a) / 0.2e1) * a + sin(pi * j * (-DeltaT * sigma ^ 2 + 0.2e1 * DeltaT * r - 0.4e1 * a + 0.4e1 * l) / (-b + a) / 0.2e1) * b + sin(pi * j * (-DeltaT * sigma ^ 2 + 0.2e1 * DeltaT * r - 0.4e1 * a + 0.4e1 * h) / (-b + a) / 0.2e1) * a - sin(pi * j * (-DeltaT * sigma ^ 2 + 0.2e1 * DeltaT * r - 0.4e1 * a + 0.4e1 * h) / (-b + a) / 0.2e1) * b) / pi / j / 0.4e1);
%   end
%   
%   A3kj_Matrix( 1, 1)  =  2/(b-a) * (h-l ) ;
%   A3kj_Matrix( :, 1) =   A3kj_Matrix(  :, 1)  / 2 ;  
    
%   Akj_Matrix = A3kj_Matrix;
 

%% %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%
end


