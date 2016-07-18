function data = preprocess_data(individual_data_file, choice_data_file, choice_characters_file, num_choices_per_person)
%
%data = preprocess_data(individual_data_file, choice_data_file, choice_characters_file)
%
%Preprocesses the data on individuals, their choices and the choice
%characteristics
%
%File properties:
%
%individual_data_file: a tabulated list of individual ID, individual_location_x, individual_location_y, individual
%character 1, individual character 2 ... individual character n
%
%choice_data_file: a tabulated list of individual ID, choice ID
%
%choice_characters_file: a tabulated list of choice ID, choice_location_x, choice_location_y, number of units (e.g. # houeses or schools) set to 1 if not needed), choice character 1,
%choice character 2, choice character 3... choice chracter n
%
%num_choices_per_person: number of choices to sample for each person.
%Actual selected choice is always sampled.

individual_data = load(individual_data_file);
choice_data = load(choice_data_file);
choice_characters = load(choice_characters_file);

%find individuals for whom we have choices
[individual_ids, IA, IB] = intersect(individual_data(:, 1), choice_data(:, 1));

individual_data = individual_data(IA, :);
choice_data = choice_data(IB, :);

choice_ids = choice_characters(:, 1);

data = NaN*ones(numel(individual_ids), num_choices_per_person, size(individual_data,2)+size(choice_characters,2) - 2 - 4 +1); %-2 for ids, -4 for locations, +1 for distance

for i = 1:numel(individual_ids)
    
    individual = individual_data(i, :);
    choice_made = choice_data(i, 2);
    
    %Choose random choices for comparison
    choices_to_consider = randsample(choice_ids, num_choices_per_person);
    
    %if actual choice is not included, include it 1st, else find it and
    %move to top
    idx = find(choices_to_consider==choice_made);
    if isempty(idx)
        choices_to_consider(1) = choice_made;
    else
        choices_to_consider([1, idx]) = choices_to_consider([idx, 1]);
    end
    
    
    
    choices_considered_characters = choice_characters(choices_to_consider, :);
    distance = sqrt(sum(bsxfun(@minus, choices_considered_characters(:, 3:4), individual(2:3)).^2, 2));
    log_num_units = log(choices_considered_characters(:, 2));
 
    data(i, :, :) = [log_num_units, repmat(individual(4:end), num_choices_per_person, 1), distance, choices_considered_characters(:, 5:end)]; 
    
end

    







