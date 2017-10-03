function [chi, psi] = coeff(k, c, d, a, b)
diff1 = b-a;
diff2 = d-a;
diff3 = c-a;

chi = (exp(d)-exp(c));  %holds for k = 0
psi = (d-c);            %" "

auxVar = double(diff1 ./ (k .* pi) ...
        .* ( sin(k .* pi .* diff2 ./ diff1) - ...
        sin(k .* pi .* diff3 ./ diff1) ));      

psi(2:end,:) = auxVar(2:end,:);
    
chi1 = 1 ./ ( 1 + ( k .* pi ./ diff1 ) .^ 2 );
chi2 = exp(d) .* cos(k .* pi .* diff2 ./ diff1) ...
        - exp(c) .* cos(k .* pi .* diff3 ./ diff1);
chi3 = k .* pi ./ (b - a) .* ... 
        ( exp(d) .* sin(k .* pi .* diff2 ./ diff1) ...
        - exp(c) .* sin(k .* pi .* diff3 ./ diff1) );

auxVar = chi1 .* (chi2 + chi3);

chi(2:end,:) = auxVar(2:end,:);

end