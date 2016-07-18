function lP = logmvnpdf_cholesky(X, R)
%lP = logmvnpdf_cholesky(X,R)
%
%calculates the log-probability density at X given cholesky
%decomposition of the covariance matrix R for a multivariate normal
%distribution. It is assumed that distribution is 0 mean.
  
%Richard Mann (2009)

%Options for linsolve
upper.UT = true;
lower.UT = true;
lower.TRANSA = true;

X = X(:);

lP = -0.5*(X'*linsolve(R, linsolve(R, X, lower), upper) + 2*sum(log(diag(R))) + length(R)*log(2*pi));