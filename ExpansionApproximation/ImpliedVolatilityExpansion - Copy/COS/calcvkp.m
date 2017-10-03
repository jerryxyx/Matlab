function y = calcvkp(k, b, a, strike)

[chi, psi] = coeff(k,a,0,a,b);
y = (psi - chi)*diag(strike);

end

