
function [  f_k ] = f_k_expansion(   a, b, l, h,  K,  S0,   N_FouriersTerms    )

 
%% the code computes the f_k =  

 
%%
N = N_FouriersTerms;
 
vec_k = [0:N];

f_k= zeros(N,1);


strike = log( K./ S0 );


 

for k_index=1:length(vec_k)
    k = vec_k(k_index);
    f_k( k_index )  =  real( (-b + a) * (exp(-l * ((1i) * pi * k + b - a) / (-b + a)) * k * pi - exp(-(strike * b - strike * a + (1i) * k * pi * l) / (-b + a)) * pi * k + (i) * exp(-(strike * b - strike * a + (i) * k * pi * l) / (-b + a)) * b + (-1*i) * exp(-(strike * b - strike * a + (i) * k * pi * l) / (-b + a)) * a - exp(-h * ((i) * pi * k + b - a) / (-b + a)) * k * pi + exp(-(strike * b - strike * a + (i) * k * pi * h) / (-b + a)) * pi * k + (-1*i) * exp(-(strike * b - strike * a + (i) * k * pi * h) / (-b + a)) * b + (i) * exp(-(strike * b - strike * a + (i) * k * pi * h) / (-b + a)) * a) * exp((i) * k * pi / (-b + a) * a) / ((i) * pi * k + b - a) / k / pi ) ;
end
 
    f_k( 1 )  =  real( -exp(l) + exp(strike) * l + exp(h) - exp(strike) * h) ;
 

f_k = f_k * S0  ; 
 



%  
% f2_k = zeros(N ,1);
% 
% for k_index=1:length(vec_k)
%     k = vec_k(k_index);
%     f= @(z) ( S0 .* (exp(z) - exp( strike )) .* cos(k .* pi .* (z - a) ./ (b - a)) );
%      f2_k( k_index ) =  integral(  f, l, h,  'ArrayValued',false ) ;   
% %     f_k(i) =  quad(  f, delta , b) ;  
% end
% 
%  

 

end





