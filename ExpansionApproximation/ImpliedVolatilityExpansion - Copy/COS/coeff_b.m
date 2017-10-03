function [chi, psi] = coeff_b(k, x1, x2, a, b)
% compute chi and psi for cosine method

arg2 = k .* pi * diag((x2 - a) ./ (b - a));     % arg trig func
arg1 = k .* pi * diag((x1 - a) ./ (b - a));     % arg trig func

term1 = cos( arg2 ) * diag(exp(x2));
term2 = cos( arg1 ) * diag(exp(x1));

term3 = pi * k .* sin( arg2 ) * diag(exp(x2)./ (b-a));
term4 = pi * k .* sin( arg1 ) * diag(exp(x1)./ (b-a));

chi = 1 ./ ( 1 + ((k .* pi) *diag(1./ (b - a))).^2 ) ...
    .* ( term1 - term2 + term3 - term4 );   % modify init

chi(1,:) = (exp(x2)-exp(x1)); % init chi

psi = ((sin(arg2) - sin(arg1)) ./ (k .* pi)) *diag(b-a);    % modify init

psi(1,:) = (x2-x1);           % init psi

end
