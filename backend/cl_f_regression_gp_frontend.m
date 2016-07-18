function [f, fval, W, nrawLL] = cl_f_regression_gp_frontend(log_num_units, selection_probs, f,N, ind, R)


nLL = @(params) cl_regression_f_LL_gp_frontend(log_num_units, selection_probs, params, N, ind, R);
LLforW = @(params) cl_regression_f_LL_gp_frontend(log_num_units, selection_probs, params, N, ind);

start = vertcat(f{:});


options.numDiff = 0;
options.Display = 'iter';
options.DerivativeCheck = 'off';
options.MaxIter = 200;
options.UseMex = 1;


[BETA, fval] = minFunc(nLL, start, options);

for i = 1:numel(f)
    if i == 1
        startidx = 1;
    else
        startidx = sum(N(1:i-1))+1;
    end
    f{i} = BETA(startidx:startidx+N(i)-1);
end


if nargout >= 3
    disp('Numerically estimating Hessian')
    
   [~,~,W] = autoHess(BETA,1, LLforW);
 
  
end

if nargout >= 4
    nrawLL = LLforW(BETA);
end

