function [hyp, fval] = cl_hyper_regression_gp_frontend(f,Z, N, n,hyp, Rw, covariance_function)


nLL = @(params) cl_regression_hyper_LL_gp_frontend(f, Z, N, n, params, Rw, covariance_function);


options.numDiff = 0;
options.Display = 'iter';
options.DerivativeCheck = 'off';
options.MaxIter = 50;
options.UseMex = 1;
options.Method = 'lbfgs';

if strcmp(covariance_function, 'matern')
for i = 1:numel(hyp)
    hyp{i}(end) = log(std(f{i})); %move the starting output scale to somewhere reasonable - only applies when output scale exists!
end
end

start = vertcat(hyp{:});


[hypout, fval] = minFunc(nLL, start, options);

for i = 1:numel(hyp)
if i == 1
        startidx = 1;
else
        startidx = cumsum(n(1:i-1))+1;
end
hyp{i} = hypout(startidx:startidx+n(i)-1);
end





