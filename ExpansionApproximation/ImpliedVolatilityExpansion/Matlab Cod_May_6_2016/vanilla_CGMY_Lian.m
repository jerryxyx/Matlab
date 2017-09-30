% vanilla_CGMY
% 
% 
% 
% 
% 
% 




S = 100; 
K  = 100;
r =0.05;
d = 0.02;
T = 1;
 


%model parameters CGMY
C = 4; G = 50; M = 60; Y = 0.7;
modelCGMY.ID = 'CGMY';          %model identifyer
modelCGMY.params = [ C;        %C 
                     G;          %G
                     M;          %M
                     Y];         %Y
              
% modelCGMY.params = [1.0;        %C 
%                   5.0;          %G
%                   5.0;          %M
%                   1.5];         %Y
%               
              

[price_tmp, delta_tmp, gamma_tmp] = CosineMethodCallPricingFFT( 500, 20 , modelCGMY ,S,K,T,r,d) ;



modelNIG.ID = 'NIG';          %model identifyer
    params = [ 15, -5, 0.5];
modelNIG.params = params;

[price_tmp, delta_tmp, gamma_tmp] = CosineMethodCallPricingFFT( 500, 20 , modelNIG ,S,K,T,r,d) ;







%model parameters Heston
modelHeston.ID = 'Heston';      %model identifyer
modelHeston.params = [0.04;     %vInst 
                      0.04;     %vLong
                       1.5;     %kappa
                       0.5;     %omega
                     -0.8];     %rho             

[price_tmp, delta_tmp, gamma_tmp] = CosineMethodCallPricingFFT( 500, 20 , modelHeston ,S,K,T,r,d) 
                 
                 
                 
                 