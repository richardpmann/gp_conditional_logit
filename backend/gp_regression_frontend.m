function [fstar, sestar, fval, X,Y, hyp, neg_marginal_log_likelihood,params_record] = gp_regression_frontend(individual_data_file, choice_data_file, choice_characters_file, num_choices_per_person, em_max)
%
%gp_regression_frontend(individual_data_file, choice_data_file, choice_characters_file, num_choices_per_person)
%
%
%Frontend for the GP conditional logit model. Takes text files specifying
%individual characteristics, choices made and choice characteristics
%
%
%File properties:
%
%individual_data_file: a tabulated list of individual ID, individual_location_x, individual_location_y, individual
%character 1, individual character 2 ... individual character n
%
%choice_data_file: a tabulated list of individual ID, choice ID
%
%choice_characters_file: a tabulated list of choice ID, redundancy (number of identical choices in this unit such as number of houses in neighbourhood), choice_location_x, choice_location_y, choice character 1,
%choice character 2, choice character 3... choice character n
%
%num_choices_per_person: number of choices to sample for each person.
%Actual selected choice is always sampled.
%
%em_max: the number of Expectation Maximisation iterations to perform
%(default = 25)

%covariance_function = 'linear';
%covariance_function = 'matern';
covariance_function = 'sum';
if ~exist('em_max') || isempty(em_max)
    em_max = 25;
    
end
% if strcmp(covariance_function, 'linear')
%         em_max=1; %the linear covariance has no hyperparameters so no need to EM
% end
%Options for linsolve
upper.UT = true;
lower.UT = true;
lower.TRANSA = true;

%preprocess data into required format
disp('Pre-processing data')
data = preprocess_data(individual_data_file, choice_data_file, choice_characters_file, num_choices_per_person);

%load original data files to get key properties
individual_data = load(individual_data_file);
choice_characters = load(choice_characters_file);


disp('Setting up Gaussian process structure');
num_ind_chars = size(individual_data, 2)-3; %-3 for id and location columns
num_choice_chars = size(choice_characters, 2)-4+1; %-3 for id num units and location columns + 1 for distance

%Model is a sum of GPs, one for each choice char, function of the
%choice_char and all individual chars

%Define GPs. each starts with hyper-parameters set at log(1) for all
%scales. There is a scale for each ind_char and 1 for the output scale.

for i = 1:num_choice_chars
    if strcmp(covariance_function, 'matern')
         hyp{i} = [log(1*ones(num_ind_chars+1, 1)); log(1)]; %+1 for the choice characters and then output scale
    elseif strcmp(covariance_function, 'linear')
         hyp{i} = [log(1*ones(num_ind_chars+1, 1))]; %+1 for the choice characters, no output scale
    elseif strcmp(covariance_function, 'sum')
        hyp{i} = [log(1*ones(num_ind_chars+1, 1)); log(1);log(1*ones(num_ind_chars+1, 1))];
        covsumfn = @(varargin) covSum({{@covMaternard, 3}, @covLINard_nbdseg},varargin{:});
         
    else
        disp('Choose a valid covariance')
        
    end
end

%Individual factors. 
for i = 1:num_ind_chars
ind_factor{i} = data(:, :, i+1);

%whiten to zero mean, unit std
ind_whitening_factors{i} = [mean(ind_factor{i}(:)), std(ind_factor{i}(:))];
ind_factor{i} = (ind_factor{i}-ind_whitening_factors{i}(1))/ind_whitening_factors{i}(2);

%discretise to max 10 unique values
if numel(unique(ind_factor{i}(:))) > 10
roundTargets = quantile(ind_factor{i}(:), 0.05:0.1:0.95);
%roundTargets = linspace(min(ind_factor{i}(:)), max(ind_factor{i}(:)), 10);

ind_factor{i} = reshape(interp1(roundTargets,roundTargets,ind_factor{i}(:),'nearest', 'extrap'), size(ind_factor{i}));

end

end


%Choice factors

log_num_units = data(:,:, 1);
selection_probs = ((num_choices_per_person-1)/(size(choice_characters, 1)-1))*ones(size(data(:, :, 1)));
selection_probs(:, 1, :) = 1;


for i = 1:num_choice_chars
    choice_factor{i} = data(:, :, num_ind_chars+1+i);
    
    %whiten to zero mean, unit std
    choice_whitening_factors{i} = [mean(choice_factor{i}(:)), std(choice_factor{i}(:))];
    choice_factor{i} = (choice_factor{i}-choice_whitening_factors{i}(1))/choice_whitening_factors{i}(2);
    
    %discretise to max 10 unique values
    if numel(unique(choice_factor{i}(:))) > 10
    roundTargets = quantile(choice_factor{i}(:), 0.05:0.1:0.95)+eps*rand(1, 10);
    roundTargets=sort(roundTargets, 'ascend');
    %roundTargets = linspace(min(choice_factor{i}(:)), max(choice_factor{i}(:)), 10);

    choice_factor{i} = reshape(interp1(roundTargets,roundTargets,choice_factor{i}(:),'nearest', 'extrap'), size(choice_factor{i}));

    end
    
end

%find unique combinations of individual and choice chars and their indices

for i = 1:num_choice_chars
   
    buildingZ = NaN*ones(size(ind_factor{1}(:), 1), num_ind_chars+1);
    for j = 1:num_ind_chars
        buildingZ(:, j) = ind_factor{j}(:);
    end
    buildingZ(:, num_ind_chars+1) = choice_factor{i}(:);
   
   
[Z{i}, ~, ind{i}] = unique(buildingZ, 'rows');

N(i) = size(Z{i}, 1); %number of points on GP i
if strcmp(covariance_function, 'matern')
n(i) = size(Z{i}, 2)+1; %number of hyperparams needed
elseif strcmp(covariance_function, 'linear')
n(i) = size(Z{i}, 2); %number of hyperparams needed
elseif strcmp(covariance_function, 'sum')
n(i) = size(Z{i}, 2)*2 + 1; %number of hyperparams needed
end
f{i} = zeros(N(i), 1); %(rand(N(i), 1)-0.5)*0.1; %initiate the mean of GP i

end



tic
disp(['There are ', num2str(num_ind_chars), ' individual characteristics and ', num2str(num_choice_chars), ' choice characteristics']);
disp('Performing GP conditional logit regression');

fval_for_convergence = inf;
for em_loop = 1:em_max
disp(['EM iteration: ', num2str(em_loop)])

%covariances

for i = 1:num_choice_chars
    if strcmp(covariance_function, 'matern')
        K{i} = covMaternard(3, hyp{i}, Z{i});
        R{i} = chol(K{i});
    elseif strcmp(covariance_function, 'linear')
        K{i} = covLINard_nbdseg(hyp{i}, Z{i});
        K{i} = K{i} + eye(size(K{i}))*0.001;
        R{i} = chol(K{i});
        %R{i} = bsxfun(@times, [Z{i}(:, 1:(size(Z{i}, 2)-1)), ones(size(Z{i}(:, 1)))], Z{i}(:, end))';
    elseif strcmp(covariance_function, 'sum')
       
        K{i} = covsumfn(hyp{i}, Z{i});
         K{i} = K{i} + eye(size(K{i}))*0.001;
        R{i} = chol(K{i});
    else
        disp('Choose a valid covariance function');
    end
  

end



disp('Learning utility function')
[fout, fval(em_loop), Wall, nrawLL] = cl_f_regression_gp_frontend(log_num_units, selection_probs, f,N, ind, R);





Wall = Wall + 0.001*eye(size(Wall)); %jitter to deal with numerical precision.
Rw = chol(Wall);


disp('Learning automatic relevance detection parameters')
[hypout, hyperfval(em_loop)] = cl_hyper_regression_gp_frontend(fout, Z, N, n, hyp, Rw, covariance_function);


neg_marginal_log_likelihood(em_loop) = (hyperfval(em_loop)+nrawLL);
%if the loglikelihood doesn't improve by at least 1, we declare convergence
%and break the EM loop
if fval_for_convergence-neg_marginal_log_likelihood(em_loop) > 1
    fval_for_convergence = neg_marginal_log_likelihood(em_loop);
    f = fout;
    hyp = hypout;
else
    break;
end

params_record(:, em_loop) = vertcat(hyp{:});

toc
end
    

% for i = 1:num_choice_chars
% K{i} = covMaternard(3, hyp{i}, Z{i});
% iK{i} = inv(K{i});
% end

%S = inv(Wall+blkdiag(iK{:}));

Expansion_Factor = 20; %Expands results to a denser grid for plotting (>=50 known to break on 16GB RAM)

for i = 1:num_ind_chars
    X{i} = unique(ind_factor{i}(:));
    if numel(X{i}) == 10
        X{i} = linspace(min(ind_factor{i}(:)), max(ind_factor{i}(:)), Expansion_Factor)';
    end
end

for i = 1:num_choice_chars
    Y{i} = unique(choice_factor{i}(:));
    if numel(Y{i}) == 10
        Y{i} = linspace(min(choice_factor{i}(:)), max(choice_factor{i}(:)), Expansion_Factor)';
    end
    
end

for i = 1:numel(f)

    incell = [X, Y(i)];
    outcell = cell(1,numel(incell));
    [outcell{:}] = ndgrid(incell{:});
    Zstar = [];
    Zstardim = [];
    for j = 1:numel(outcell)
        Zstar = [Zstar, outcell{j}(:)];
        Zstardim(j) = numel(incell{j});
    end
    
 if strcmp(covariance_function, 'matern')
    K = covMaternard(3, hyp{i}, Z{i});
    R = chol(K);
    Kstar = covMaternard(3, hyp{i}, Zstar, Z{i});
    Kstarstar = covMaternard(3, hyp{i}, Zstar);
 elseif strcmp(covariance_function, 'linear')
    K = covLINard_nbdseg(hyp{i}, Z{i});
    K = K + eye(size(K))*0.001;
    R = chol(K);
    %R = bsxfun(@times, [Z{i}(:, 1:(size(Z{i}, 2)-1)), ones(size(Z{i}(:, 1)))], Z{i}(:, end))';

    Kstar = covLINard_nbdseg(hyp{i}, Zstar, Z{i});
    Kstarstar = covLINard_nbdseg(hyp{i}, Zstar);
 
 elseif strcmp(covariance_function, 'sum')
    K = covsumfn(hyp{i}, Z{i});
    R = chol(K);

    Kstar = covsumfn(hyp{i}, Zstar, Z{i});
    Kstarstar = covsumfn(hyp{i}, Zstar);   
 else
        disp('Choose a valid covariance function');
 end
   
    

    

    fstar{i} = Kstar*linsolve(R, linsolve(R, f{i}, lower), upper);
    sestar{i} = sqrt(diag(Kstarstar - Kstar*linsolve(R, linsolve(R, Kstar', lower), upper)));

    fstar{i} = reshape(fstar{i}, Zstardim);
    sestar{i} = reshape(sestar{i}, Zstardim);
    
end

for i = 1:num_ind_chars
    X{i} = X{i}*ind_whitening_factors{i}(2) + ind_whitening_factors{i}(1);
end

for i = 1:num_choice_chars
    Y{i} = Y{i}*choice_whitening_factors{i}(2) + choice_whitening_factors{i}(1);
end


 




