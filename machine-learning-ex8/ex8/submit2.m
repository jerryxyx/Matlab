n_u = 3; n_m = 4; n = 5;
  X = reshape(sin(1:n_m*n), n_m, n);
  Theta = reshape(cos(1:n_u*n), n_u, n);
  Y = reshape(sin(1:2:2*n_m*n_u), n_m, n_u);
  R = Y > 0.5;
  pval = [abs(Y(:)) ; 0.001; 1];
  Y = (Y .* double(R));  % set 'Y' values to 0 for movies not reviewed
  yval = [R(:) ; 1; 0];
  params = [X(:); Theta(:)];

    [mu sigma2] = estimateGaussian(X);
    out = sprintf('%0.5f ', [mu(:); sigma2(:)]);