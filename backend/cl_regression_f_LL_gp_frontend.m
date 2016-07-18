function [nLL, dnLL] = cl_regression_f_LL_gp_frontend(log_num_units, selection_probs, fin,N, ind, R)

%Options for linsolve
upper.UT = true;
lower.UT = true;
lower.TRANSA = true;

%Parcel input f into cell
for i =1:numel(N)
    if i == 1
        startidx = 1;
    else
        startidx = sum(N(1:i-1))+1;
    end
    f{i} = fin(startidx:startidx+N(i)-1);
end

LL = 0;
dLL = zeros(sum(N), 1);

if exist('R', 'var') && ~isempty(R)
%Prior on latent function

for i = 1:numel(N)
priorLL = logmvnpdf_cholesky(f{i}, R{i}); %this does the regularisation using the GP prior
LL =  LL + priorLL;

if i == 1
        startidx = 1;
else
        startidx = sum(N(1:i-1))+1;
end
    
dLL(startidx:startidx+N(i)-1) = -linsolve(R{i}, linsolve(R{i}, f{i}, lower), upper);

end

end

Nind = size(log_num_units, 1);
Nchoice = size(log_num_units, 2);


U = zeros(Nind, Nchoice);

%Correction for variable size nbd
U = U + log_num_units;

for i = 1:numel(f)

U = U + reshape(f{i}(ind{i}), size(U));

end

U = U - log(selection_probs); % reweighting for selections

P = exp(bsxfun(@minus, U, max(U, [], 2)));
P = bsxfun(@times, P, 1./sum(P, 2));

rawLL = sum(log(P(:,1))); %this is the log probability of the data, given f
LL = LL + rawLL;

%Derivatives

dlPkdUj = -P;
dlPkdUj(:, 1) = dlPkdUj(:, 1) + 1;

for i = 1:numel(f)
    
    if i == 1
        startidx = 1;
    else
        startidx = sum(N(1:i-1))+1;
    end
   
    dLLi = accumarray(ind{i}, dlPkdUj(:), [N(i), 1]);
   
    dLL(startidx:startidx+N(i)-1) = dLL(startidx:startidx+N(i)-1) + dLLi;
    
end

nLL = -LL;
dnLL = -dLL;





    
    
    
    
    
    
    











