

function [ Akj_Matrix ] = A_Maxtrix_LevyProcesses_version2(  a, b, l,  h,   DeltaT,  N_FouriersTerms,   model, OptionParameters, Akj_Matrix_SecondPart   )

%% The code generates the A_jkn Matrix


%% xa, xb are the boundaries of the Fourier expansion to truncate the integration;
%% c and d are the boundaries computing the integration A_j,k,n


%%
 


%% Declair the total number of terms of the Fourier expansions
%%
N = N_FouriersTerms;


vec_j = [0:N];
vec_k = [0:N];


Akj_Matrix= zeros(N+1, N+1 );




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
    Akj_Matrix(:, k_index) =  real(  exp_psi_fun_vector( k_index )  * Akj_Matrix_SecondPart(:, k_index)   );
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


% COS_k_vector_l = zeros(1, length(vec_k));
% SIN_k_vector_l = zeros(1, length(vec_k));
% COS_k_vector_l = cos(pi * vec_k * (- l + a) / (-b + a));
% SIN_k_vector_l  = sin(pi * vec_k * (- l + a) / (-b + a));
% 
% COS_k_vector_u = zeros(1, length(vec_k));
% SIN_k_vector_u = zeros(1, length(vec_k));
% COS_k_vector_u = cos(pi * vec_k * (- h + a) / (-b + a));
% SIN_k_vector_u = sin(pi * vec_k * (- h + a) / (-b + a));
% 
% EXP_j_vector_u = zeros(1, length(vec_j));
% EXP_j_vector_u = exp( 1i * pi * vec_j * (- h + a) / (-b + a));
%  
% EXP_j_vector_l = zeros(1, length(vec_j));
% EXP_j_vector_l = exp( 1i * pi * vec_j * (- l + a) / (-b + a));
% 
% 
% for k_index=1:length(vec_k)
%     k = vec_k(k_index);
%      for j_index=1:length(vec_j)
%         j = vec_j(j_index);   
%         if( j~=k)
% %             Akj_Matrix(k_index, j_index)  = 2/(b-a) * real(  exp( DeltaT * psi_fun_vector(j_index)  ) * (-((i) * j * exp((i) * pi * j * (-l + a) / (-b + a)) * cos(pi * k * (-l + a) / (-b + a)) + k * exp((i) * pi * j * (-l + a) / (-b + a)) * sin(pi * k * (-l + a) / (-b + a)) + (-1*i) * j * exp((i) * pi * j * (-h + a) / (-b + a)) * cos(pi * k * (-h + a) / (-b + a)) - k * exp((i) * pi * j * (-h + a) / (-b + a)) * sin(pi * k * (-h + a) / (-b + a))) * (-b + a) / pi / (j ^ 2 - k ^ 2) ) );
% %            Akj_Matrix(k_index, j_index)  = 2 / pi / (j ^ 2 - k ^ 2)  * real( exp_psi_fun_vector(j_index) * (  (1i*j* COS_k_vector_l(k_index) + k*SIN_k_vector_l(k_index) )* EXP_j_vector_l(j_index)  - ( 1i*j* COS_k_vector_u(k_index) + k*SIN_k_vector_u(k_index)  )* EXP_j_vector_u(j_index) )   ) ;  
%             Akj_Matrix(k_index, j_index)  =   real( exp_psi_fun_vector(j_index) * Akj_Matrix_SecondPart(k_index, j_index)  );
%         end
%         
%     end
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%  Akj_Matrix(1, 1) =  2/(b-a) * real(  exp( DeltaT * psi_fun( k*pi/(b-a) ) ) * (h-l ) );
  Akj_Matrix(1, 1) =    2/(b-a) * ( h-l ) ;
 %% Because the A_j is the coefficient of the Fourier series, A_0 should be divided by 2; 
 Akj_Matrix( :, 1) =  Akj_Matrix( :, 1)  / 2 ;  
 
 

%% %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%
end


