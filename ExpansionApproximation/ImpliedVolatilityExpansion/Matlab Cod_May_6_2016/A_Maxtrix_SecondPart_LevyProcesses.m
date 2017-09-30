

function [ Akj_Matrix_SecondPart ] = A_Maxtrix_SecondPart_LevyProcesses(  a, b, l,  h,  N_FouriersTerms   )

%% The code generates the A_jkn Matrix


%% xa, xb are the boundaries of the Fourier expansion to truncate the integration;
%% c and d are the boundaries computing the integration A_j,k,n


%%
 


%% Declair the total number of terms of the Fourier expansions
%%
N = N_FouriersTerms;


vec_j = [0:N];
vec_k = [0:N];


Akj_Matrix_SecondPart= zeros(N+1,N+1);




%% ksi( ) is a function in the characteristic function exponential. 
%% If k == j then ...
for k_index=1:length(vec_k)
    k = vec_k(k_index);
    Akj_Matrix_SecondPart(k_index, k_index) =  2/(b-a) * (   ( (((i) * a + (-1*i) * b) * sin(pi * k * (-l + a) / (-b + a)) ^ 2 + cos(pi * k * (-l + a) / (-b + a)) * (-b + a) * sin(pi * k * (-l + a) / (-b + a)) + ((-1*i) * a + (i) * b) * sin(pi * k * (-h + a) / (-b + a)) ^ 2 - cos(pi * k * (-h + a) / (-b + a)) * (-b + a) * sin(pi * k * (-h + a) / (-b + a)) - pi * k * (l - h)) / pi / k / 2 )   );
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
            Akj_Matrix_SecondPart(k_index, j_index)  = 2 / pi / (j ^ 2 - k ^ 2)  * ( (  (1i*j* COS_k_vector_l(k_index) + k*SIN_k_vector_l(k_index) )* EXP_j_vector_l(j_index)  - ( 1i*j* COS_k_vector_u(k_index) + k*SIN_k_vector_u(k_index)  )* EXP_j_vector_u(j_index) )   ) ;  
        end
        
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 Akj_Matrix_SecondPart(1, 1) =    2/(b-a) * ( h-l ) ;



end


