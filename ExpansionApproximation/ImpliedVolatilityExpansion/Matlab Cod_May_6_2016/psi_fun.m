

function [ f  ]  = psi_fun( lambda, model, OptionParameters  )
%% The function ksi in the function B7.



% S0 = OptionParameters(1);
% K = OptionParameters(2);
% T = OptionParameters(3);
r = OptionParameters(4);
d = OptionParameters(5);

% CallPutIndicator    = OptionParameters(6);
% InOutIndicator      = OptionParameters(7);
% LowerBoundary       = OptionParameters(8);
% UpBarrier           = OptionParameters(9);
% NumObservations     = OptionParameters(10);





% mu =  ( r -div ) - sigma^2 /2  - integration; 
 
%  f = -  lambda^2 * sigma^2 /2  + i * lambda * mu +  integration;
 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diffusion component
if strcmp(model.ID,'BlackScholes')
    sigma = model.params(1);
 mu  =  ( r - d ) - sigma^2 /2 ; 
 f = -  lambda^2 * sigma^2 /2  + i * lambda * mu ;

 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Black-Scholes  %%     
    %% Note the parameter: 
%         u = lambda;  T =1;  
% %        sigma = params(1); a = params(2); b = params(3); lambda = params(4);
%       f  = Black_characteristicFn(u,T,r,d, model.params);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 %% The NIG Model

%  f= delta *( sqrt( alpha^2 - (beta+i*lambda)^2 ) - sqrt(alpha^2 - beta^2 )  )




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CGMY
        elseif strcmp(model.ID,'CGMY')
           
% C = 4; G = 50; M = 60; Y = 0.7;
% C =  1.0;        %C 
% G =  5.0;        %G
% M =  5.0;        %M
% Y =  1.5;         %Y
              
C =  model.params(1);        
G =  model.params(2);         
M =  model.params(3);        
Y =  model.params(4);         
T = 1;
    omega = -C*gamma(-Y)*((M-1)^Y-M^Y+(G+1)^Y-G^Y);
    %tmp = exp(C*T*gamma(-Y)*((M-i*u).^Y-M^Y+(G+i*u).^Y-G^Y));
    %phi = exp( i*u*(lnS + (r-d+m)*T)).*tmp;
    tmp = C*T*gamma(-Y)*((M-1i* lambda ).^Y-M^Y+(G+1i* lambda ).^Y-G^Y);
    f = 1i* lambda *((r-d+omega)*T) + tmp;
%   f =  CGMY_characteristicFn(lambda,T,r,d,model.params)  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VG Processes  %%

%% Below is the same the code of the function phi = VG_characteristicFn(u,T,r,d,params)

%     sigma = params(1); nu = params(2); theta = params(3);    
%     u = lambda;  T =1; d = div;
%     
%     params = [4, 12, 18];
%     
%     c = params(1); g = params(2); m = params(3);
%     nu = 1/c;
%     theta = (1/m-1/g)/nu;
%     sigma = sqrt(((1/g+theta*nu/2)^2-theta^2*nu^2/4)*2/nu);
%     
%     omega = (1/nu)*( log(1-theta*nu-sigma*sigma*nu/2) );
%     tmp = 1 - 1i * theta * nu * u + 0.5 * sigma * sigma * u .* u * nu;
%     %tmp = tmp.^(-T / nu);
%     %phi = exp( i * u * (lnS + (r + omega - d) * T )) .* tmp;
%     %phi = 1i * u * (lnS + (r + omega - d) * T ) - T*log(tmp)/nu;
%     phi = 1i * u * (r-d+omega)* T  - T*log(tmp)/nu;
% 
%     f = phi ;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    

 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NIG Processes  %%   
%                 u = lambda;  T =1; d = div;
%                 params = [ 15, -5, 0.5];
%                 alfa = params(1); beta = params(2); delta = params(3);
% 
%                 omega = delta*(sqrt(alfa*alfa-(beta+1)^2)-sqrt(alfa*alfa-beta*beta));
%                 %tmp = exp(i*u*mu*T-delta*T*(sqrt(alfa*alfa-(beta+i*u).^2)-sqrt(alfa*alfa-beta*beta)));
%                 %phi = exp( i*u*(lnS + (r-d+m)*T)).*tmp;
%                 tmp = -delta*T*(sqrt(alfa*alfa-(beta+1i*u).^2)-sqrt(alfa*alfa-beta*beta));
%                 phi = 1i*u*((r-d+omega)*T) + tmp;
%                 f = phi ;   
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Merton Jump + diffusion Processes  %%  
        elseif strcmp(model.ID,'Merton')
        u = lambda;  T =1;  
        
    %% Note the parameter: exp(m+.5*b^2)-1.0 = a 
%         params = [ 0.1, -0.05, 0.086, 3];
%         sigma = params(1); a = params(2); b = params(3); lambda = params(4);
%         model.params(2) = log(1+ model.params(2))- 0.5 * model.params(3)^2 ;
        phi = Merton_characteristicFn(u,T,r,d, model.params);
        f = phi ; 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
    

    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Kou Jump model  %%     
    %% Note the parameter: 
        

        elseif strcmp(model.ID,'Kou')
%         params = [ 0.1, 3, 0.3, 40, 12 ];  %% the parameters in Feng Math Finance 2007
%       params = [ 0.353, 6.37, 0.6, 10, 5.712 ];  %% the parameters in Banking Finance 2013  

%        params = [ 0.212, 2.29, 0.6, 10, 5.712];  %% Table 4 parameters  the parameters in Banking Finance 2013  
%        params = [ 0.0, 2.29, 0.6, 10, 5.712];  %% Table 4 parameters  the parameters in Banking Finance 2013  
%     params = [ 0.3, 1, 0.5, 30, 30 ];  %% Table 1 parameters  the parameters in Banking Finance 2013  
%      params = [  1.011512676508229, 12, 0.5, 10, 10 ];  %% Table 2 parameters  the parameters in Banking Finance 2013  
%     sigma  = params(1);  lambda = params(2);  pUp   = params(3);  %% p
%     mUp   = params(4); %% eta_1
%     mDown = params(5); %% eta_3        
    u = lambda;  T = 1;  
    phi = Kou_characteristicFn(u,T,r,d, model.params);
    f   = phi ;
end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     phi = Kou_characteristicFn(u,T,r,d,params)
    
    
end



% Black Scholes    
function phi = Black_characteristicFn(u,T,r,d,params)
    sigma = params(1);
    %phi = exp(i*u*(lnS+(r-d-0.5*sigma*sigma)*T) - 0.5*sigma*sigma*u.*u*T);
     %phi = 1i*u*(lnS+(r-d-0.5*sigma*sigma)*T) - 0.5*sigma*sigma*u.*u*T;
      omega = -0.5*sigma*sigma;
      phi = 1i*u*((r-d+omega)*T) - 0.5*sigma*sigma*u.*u*T;
end

% Merton Jump Diffusion
function phi = Merton_characteristicFn(u,T,r,d, params)
    sigma = params(1); a = params(2); b = params(3); lambda = params(4);
    phi = Black_characteristicFn(u,T,r,d,sigma) + LogNormalJump_characteristicFn(u,T,[a;b;lambda]);
end

% LogNormalJump for Merton and Bates
function phiJump = LogNormalJump_characteristicFn(u,T,params)

    a = params(1); b = params(2); lambda = params(3);

    %phiJump = exp(lambda*T*(-a*u*i + (exp(u*i*log(1.0+a)+0.5*b*b*u*i.*(u*i-1.0))-1.0)));
    %phiJump = lambda*T*(-a*u*1i + (exp(u*1i*log(1.0+a)+0.5*b*b*u*1i.*(u*1i-1.0))-1.0));
    
    phiJump = lambda*T*( exp(u*1i* a - 0.5*b*b*u*u ) -1.0 ) -  1i*u*( lambda*T*( exp(1* a + 0.5*b*b ) -1.0 ) );
end



% Kou Jump Model
function phi = Kou_characteristicFn(u,T,r,d,params)

    sigma  = params(1);  
    lambda = params(2);  
    pUp   = params(3);  %% p
    mUp   = params(4); %% eta_1
    mDown = params(5); %% eta_3
    
    comp            = @(x) pUp*mUp./(mUp-x) + (1-pUp)*mDown./(mDown+x)-1;
    phi  = Black_characteristicFn(u,T,r,d,sigma) + lambda*comp(u*i)  -  u*i*lambda*comp(1);
end



% CGMY
function phi = CGMY_characteristicFn(u,T,r,d,params)

    C = params(1); G = params(2); M = params(3); Y = params(4);

    omega = -C*gamma(-Y)*((M-1)^Y-M^Y+(G+1)^Y-G^Y);
    %tmp = exp(C*T*gamma(-Y)*((M-i*u).^Y-M^Y+(G+i*u).^Y-G^Y));
    %phi = exp( i*u*(lnS + (r-d+m)*T)).*tmp;
    tmp = C*T*gamma(-Y)*((M-1i*u).^Y-M^Y+(G+1i*u).^Y-G^Y);
    phi = 1i*u*((r-d+omega)*T) + tmp;
end



% Variance Gamma
function phi = VG_characteristicFn(u,T,r,d,params)

    %sigma = params(1); nu = params(2); theta = params(3);
    c = params(1); g = params(2); m = params(3);
    nu = 1/c;
    theta = (1/m-1/g)/nu;
    sigma = sqrt(((1/g+theta*nu/2)^2-theta^2*nu^2/4)*2/nu);
    
    omega = (1/nu)*( log(1-theta*nu-sigma*sigma*nu/2) );
    tmp = 1 - 1i * theta * nu * u + 0.5 * sigma * sigma * u .* u * nu;
    %tmp = tmp.^(-T / nu);
    %phi = exp( i * u * (lnS + (r + omega - d) * T )) .* tmp;
    %phi = 1i * u * (lnS + (r + omega - d) * T ) - T*log(tmp)/nu;
    phi = 1i * u * (r-d+omega)* T  - T*log(tmp)/nu;
    
end




% Normal Inverse Gaussian
function phi = NIG_characteristicFn(u,T,r,d,params)

    alfa = params(1); beta = params(2); delta = params(3);

    omega = delta*(sqrt(alfa*alfa-(beta+1)^2)-sqrt(alfa*alfa-beta*beta));
    %tmp = exp(i*u*mu*T-delta*T*(sqrt(alfa*alfa-(beta+i*u).^2)-sqrt(alfa*alfa-beta*beta)));
    %phi = exp( i*u*(lnS + (r-d+m)*T)).*tmp;
    tmp = -delta*T*(sqrt(alfa*alfa-(beta+1i*u).^2)-sqrt(alfa*alfa-beta*beta));
    phi = 1i*u*((r-d+omega)*T) + tmp;
end


