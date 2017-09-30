function price = COS_FUNC_Call_NoPCP(rnCHF,r,T,S_0,W,N,aC,bC)
%Support on [aC,bC]


x = log(S_0/W);
b=x+bC; a=aC+x; 

U = zeros(1,N);
U(1) = 1/(b-a)*(exp(b)-1-b); %Already been divided by 2
U(2:N) = 2/(b-a)*(1./(1+((1:N-1)*pi/(b-a)).^2).*(cos((1:N-1)*pi)*exp(b)-cos((1:N-1)*pi*a/(b-a)) + pi/(b-a)*(1:N-1).*sin((1:N-1)*pi*a/(b-a))) - sin((1:N-1)*pi*a/(b-a))*(b-a)./((1:N-1)*pi));

price = exp(-r*T)*W*real(sum(U.*rnCHF((0:N-1)*pi/(b-a)).*exp(1i*(0:N-1)*pi*(x-a)/(b-a))));
end
