
function [  A_j, A_Delta_j ] = A_j_Vector(   a, b, x0, DeltaT,   N_FouriersTerms ,  model, OptionParameters , bool_Comp_Greeks   )


%% the code computes the f_k =


%%
N = N_FouriersTerms;

vec_k = [0:N];

A_j  = zeros(N,1);
A_Delta_j  = zeros(N,1);  %% Used to compute the delta of a barrier option


for k_index=1:length(vec_k)
    k = vec_k(k_index);
    temp = exp( (- 1i *  k * a * pi) / (b-a) )  *  ( exp(1i *x0 * k* pi /(b-a) )  *  exp( DeltaT * psi_fun( k* pi /(b-a) ,  model, OptionParameters )  )  );
    A_j( k_index )  = 2/ (b-a)  * real( temp  );
    if bool_Comp_Greeks==true
        A_Delta_j( k_index )  = 2/ (b-a)  * real( temp * 1i* k* pi /(b-a) )  ;
    end
end
A_j( 1 )   =   1 / (b-a) ;


A_Delta_j  = A_Delta_j  * exp(-x0) ;
A_Delta_j( 1 )  =  0 ;

%
% A2_j = zeros(N ,1);
%
% A3_j = zeros(N ,1);
% [ S0, strike, r, sigma, T, l, n ] = inputParameters( );
%
%     densityNormal= @(z) (  1./ ( sqrt( 2.*pi* sigma.^2 .* DeltaT ) ) .*  (  exp(  - ( z- (x0 + ( r - sigma.^2/2 ) * DeltaT )  ).^2 ./ (2 * sigma.^2 .* DeltaT )   )  )    );
%
% for k_index=1:length(vec_k)
%     k = vec_k(k_index);
%     f= @(z) ( 2 ./ (b-a)  .* ( densityNormal(z) ) .* cos(k .* pi .* (z - a) ./ (b - a)) );
% %     A2_j( k_index ) =  integral(  f, a, b,  'ArrayValued',false ) ;
%
%     x =0;
%    A3_j( k_index ) =  2 ./ (b-a) .* exp(-k ^ 2 * pi ^ 2 * sigma ^ 2 * DeltaT / (-b + a) ^ 2 / 0.2e1) * cos((a - x - (r - sigma ^ 2 / 0.2e1) * DeltaT) / (-b + a) * pi * k);
% %     f_k(i) =  quad(  f, delta , b) ;
% end
%
%


end





