function [nLL, dnLL] = cl_regression_hyper_LL_gp_frontend(f, Z, N, n, hypin, Rw, covariance_function)

%Parcel hypin into cell

for i = 1:numel(n)
if i == 1
        startidx = 1;
else
        startidx = cumsum(n(1:i-1))+1;
end
hyp{i} = hypin(startidx:startidx+n(i)-1);
end


%Options for linsolve
upper.UT = true;
lower.UT = true;
lower.TRANSA = true;

trace_product = @(A, B)  A(:)'*B(:);

dLL =zeros(1, numel(hypin));

%covariances
for i = 1:numel(Z)
    if strcmp(covariance_function, 'matern')
        K{i} = covMaternard(3, hyp{i}, Z{i});
        K{i} = K{i} + eye(size(K{i}))*0.001;
      
        R{i} = chol(K{i});
        
    elseif strcmp(covariance_function, 'linear') 
        K{i} = covLINard_nbdseg(hyp{i}, Z{i});
        K{i} = K{i} + eye(size(K{i}))*0.001;
        R{i} = chol(K{i});
          
    elseif strcmp(covariance_function, 'sum') 
         covsumfn = @(varargin) covSum({{@covMaternard, 3}, @covLINard_nbdseg},varargin{:});
        K{i} = covsumfn(hyp{i}, Z{i});
        K{i} = K{i} + eye(size(K{i}))*0.001;
        R{i} = chol(K{i});
    end
end




B = eye(size(Rw)) + Rw*blkdiag(K{:})*Rw';

iM = Rw'*linsolve(B, Rw);

LL = 0;
for i = 1:numel(n)
    
    if i == 1
        startidxn = 1;
    else
        startidxn = cumsum(n(1:i-1))+1;
    end
    
    if i == 1
        startidxN = 1;
    else
        startidxN = cumsum(N(1:i-1))+1;
    end

alpha = linsolve(R{i}, linsolve(R{i}, f{i}, lower), upper);

for ii = 1:n(i)
    
    if strcmp(covariance_function, 'matern')
      dK = covMaternard(3, hyp{i}, Z{i}, [], ii);
    elseif strcmp(covariance_function, 'linear') 
      dK = covLINard_nbdseg(hyp{i}, Z{i}, [], ii);
   elseif strcmp(covariance_function, 'sum') 
       dK = covsumfn(hyp{i}, Z{i}, [], ii);
    end
    dLL(startidxn+ii-1) =0.5*alpha'*dK*alpha - 0.5*trace_product(iM(startidxN:startidxN+N(i)-1, startidxN:startidxN+N(i)-1),dK);
end

LL = LL -0.5*(f{i}'*alpha);
end

LL = LL -0.5*logdet(B);


nLL = -LL;
dnLL = -dLL';







