function result = cvalue(x1, x2, a, b, N, V, t, r, RNchf)

Nstrike = length(x1);

exp2 = exp( 1i .* repmat((1:N)',1,Nstrike) * diag((x2 - a) ./ (b - a)) .* pi );    % init
exp1 = exp( 1i .* repmat((1:N)',1,Nstrike) * diag((x1 - a) ./ (b - a)) .* pi );    % init

m = zeros(3*N-1, Nstrike);                                        % init base

m(N,:) = 1i * pi * (x2 - x1) ./ (b - a);
m(N+1:2*N,:) = 1 ./ repmat((1:N)',1,Nstrike) .* ( exp2 - exp1 );
m(1:N-1,:) = - conj(flipud(m(N+1:2*N-1, :)));
m(2*N+1:3*N-1,:) = ( exp2(1:N-1,:) * diag(exp2(N, :)) - exp1(1:N-1,:) ...
    * diag(exp1(N,:)) ) ./ ( repmat((N+1:2*N-1)',1,Nstrike) );

Grid_j = (0:N-1)';                                          % fix grid

% compute u values
u = feval(RNchf, pi*repmat(Grid_j,1,Nstrike)*diag(1./(b-a))) .* V;
u(1,:) = 0.5*u(1,:);

m_s = [m(N:-1:1, :); zeros(1,Nstrike); m(2*N-1:-1:N+1, :)];
u_s = [u; zeros(N, Nstrike)];
m_c = m(3*N-1:-1:N, :);


zeta = -ones(2*N, Nstrike);
zeta(2 .* (1:N)' - 1,:) = 1;

fft_u_s = fft(u_s);
xi_s = ifft((fft(m_s)) .* fft_u_s);
xi_c = ifft((fft(m_c)) .* (zeta .* fft_u_s));

result = exp(-r * t) / pi .* imag( xi_s(1:N,:) + flipud(xi_c(1:N,:)) );

    
end

