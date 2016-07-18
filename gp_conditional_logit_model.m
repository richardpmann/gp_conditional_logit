%This script allows the user to specify three data files containing the
%necessary tabulated data for running the conditional logit model in Mann
%et al., to specify the maximum number of Expectation Maximisation iterations to
%run and the number of choices to compare for each individual who moves.
%
%Once the script has run, 7 outputs are produced. f is a cell of matrice of the inferred
%utiluty functions associated with each choice characteristic. se is the standard error on those estimates. fval gives the negative log probability of f at each
%Expectation-Maxmisation iteration. X and Y are cells of vectors giving the individual and choice characters corresponding to the dimesnions of each f. hyp
%gives the Gaussian process hyperparameters for each utility function, and
%params_record gives a record of all the hyperparameters over the EM
%iterations.
%
%The data files should be formulated in the following way:
%
%individual_data_file: a tabulated list of individual ID, individual_location_x, individual_location_y, individual
%character 1, individual character 2 ... individual character n
%
%choice_data_file: a tabulated list of individual ID, choice ID
%
%choice_characters_file: a tabulated list of choice ID, choice_log_num_units, choice_location_x, choice_location_y, choice character 1,
%choice character 2, choice character 3... choice chracter n


 individual_data_file = './example_data/sim/individual_data.txt';
choice_data_file = './example_data/sim/choice_data.txt';
 choice_characters_file = './example_data/sim/choice_characters_data.txt';
 %individual_data_file = './example_data/individual_data_binary_children.txt';
 %choice_data_file = './example_data/choice_data.txt';
 %choice_characters_file = './example_data/choice_characters_data.txt';

em_max =20;
num_choices_per_person = 100;

[f, se, fval, X,Y, hyp, neg_marginal_log_likelihood, params_record] = gp_regression_frontend(individual_data_file, choice_data_file, choice_characters_file, num_choices_per_person, em_max);

%save results in a temporary file: copy elsewhere to store
save latest_gp_conditional_logit_results_temp f se fval X Y hyp neg_marginal_log_likelihood params_record em_max num_choices_per_person

%%Example of plotting the results.

%Here we plot the utility of nbd
%%ethnicity v individual income, for individuals of swedish or non-swedish
%%households, for a household age 20-30 without children and a household age 30-60 with children. 
%change displaycharacter to plot other neighbourhood chracteristics
yaxistexts = {'Distance', 'Nbd % Non-Western', 'Nbd Mean Disp. Income', 'Nbd % With Children'};

displaycharacter = 1; %change this to display utility as a function of the neighbourhood characteristics above in yaxistexts

yaxistext = yaxistexts{displaycharacter};

colorscale = [min(f{displaycharacter}(:)), max(f{displaycharacter}(:))];

figure
subplot(2,2,1)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(1, :, 1,1,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(1, :, 1,1, :))' - 1.96*squeeze(se{displaycharacter}(1, :, 1,1, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(1, :, 1,1, :))' + 1.96*squeeze(se{displaycharacter}(1, :, 1,1, :)) < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');

axis xy
ylabel(yaxistext)
xlabel('Disposable income')
title('Swedish Young')
colorbar

subplot(2,2,2)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(2, :, 1,1,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(2, :, 1,1, :))' - 1.96*squeeze(se{displaycharacter}(2, :, 1,1, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(2, :, 1,1, :))' + 1.96*squeeze(se{displaycharacter}(2, :, 1,1, :))' < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');

axis xy
ylabel(yaxistext)
xlabel('Disposable income')
title('Non-Swedish Young')
colorbar

subplot(2,2,3)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(1, :, 2,2,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(1, :, 2,2, :))' - 1.96*squeeze(se{displaycharacter}(1, :, 2,2, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(1, :, 2,2, :))' + 1.96*squeeze(se{displaycharacter}(1, :, 2,2, :)) < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');
axis xy
ylabel(yaxistext)
xlabel('Disposable income')
title('Swedish Middle-Aged')
colorbar

subplot(2,2,4)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(2, :, 2,2,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(2, :, 2,2, :))' - 1.96*squeeze(se{displaycharacter}(2, :, 2,2, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(2, :, 2,2, :))' + 1.96*squeeze(se{displaycharacter}(2, :, 2,2, :)) < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');
axis xy
ylabel(yaxistext)
xlabel('Disposable income')
title('Non-Swedish Middle-Aged')
colorbar

