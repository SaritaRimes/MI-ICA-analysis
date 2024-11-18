%% Calculation of Mutual Information %%

clear all
close all
clc

addpath('FastICA_25', 'internalMethods');

mPu = 50;
snr = 2;
number_events_total = 2000000;
number_folds = 10;
number_events_fold = number_events_total/number_folds;

MI = cell(11, 2, number_folds); % column 1: no ICA, column 2: with ICA
totalMI = cell(11, 2, number_folds); % column 1: no ICA, column 2: with ICA

occupancies = [10 30 50 80];

for oc = occupancies
    fprintf(['Mean of pile-up: ' int2str(mPu) ',\t Occupancy: ' int2str(oc) '\n']);

    data_total = load(['../../../RuidoSimuladoNovoSimulador/TileCal/ruido_media' int2str(mPu) '/ruido_ocup' int2str(oc) '_' ...
                       int2str(number_events_total) 'sinais.txt']); % load noise data

    number_events = size(data_total, 1)/number_folds;
    number_samples = size(data_total, 2);
    number_intervals = 100; % number of discretization intervals
    h = (1-(-1))/number_intervals; % step of discretization
    index = (oc/10) + 1;

    for fold = 1:number_folds
        data = data_total((fold - 1)*number_events_fold + 1:fold*number_events_fold,:);        
    
        % Applying ICA to the data
        [dataICA, A, W] = fastica(data'); % apply ICA
        dataICA = dataICA'; % transpoe a matriz para ter as variaveis nas colunas
        
        % Normalizing the variables between -1 and 1
        normalized_noise = normalizeSamples(data);  % normalized variables
        normalized_noise_ica = normalizeSamples(dataICA); % normalized variables
        
        % Discretizing the variables and getting the number of samples per interval
        [discretized_noise, number_samples_interval] ...
            = discretizeSamples(normalized_noise, number_intervals, oc, fold);
        [discretized_noise_ica, number_samples_interval_ica] ...
            = discretizeSamples(normalized_noise_ica, number_intervals, oc, fold);
    
        % Calculating the marginal probabilities
        marginal_probability = number_samples_interval/number_events;
        marginal_probability_ica = number_samples_interval_ica/number_events;
            
        % Checking whether the marginal probabilities sums 1
        checkProbabilities("marginal", marginal_probability);
        checkProbabilities("marginal", marginal_probability_ica);
            
        % Calculating the joint probabilities
        joint_probability = jointProbability(discretized_noise, ...
                                             number_intervals, oc, fold);
        joint_probability_ica = jointProbability(discretized_noise_ica, ...
                                                 number_intervals, oc, fold);
    
        % Checking whether the joint probabilities sums 1
        checkProbabilities("joint", joint_probability);
        checkProbabilities("joint", joint_probability_ica);
        
        % Calculating the Mutual Information
        MI{index, 1, fold} = mutualInformation(size(discretized_noise, 2), number_intervals, ...
                                               marginal_probability, joint_probability, oc, fold);
        MI{index, 2, fold} = mutualInformation(size(discretized_noise_ica, 2), number_intervals, ...
                                               marginal_probability_ica, joint_probability_ica, oc, fold);
        
        % Calculating the crosstalk, which is the mutual information contained in all the process
        totalMI{index, 1, fold} = mutualInformationCrosstalk(MI{index, 1, fold}, number_samples);
        totalMI{index, 2, fold} = mutualInformationCrosstalk(MI{index, 2, fold}, number_samples);

        % Storing the matrices of mutual information in .txt files
        % Before ICA
        path = ['results-MI-ICA/DataNewNoise/MI_beforeICA_ocup' ...
                int2str(oc) '_it' int2str(fold) '.txt'];
        fopen(path, 'w');
        writematrix(MI{index, 1, fold}, path, 'Delimiter', 'tab');
        % After ICA
        path = ['results-MI-ICA/DataNewNoise/MI_afterICA_ocup' ...
                int2str(oc) '_it' int2str(fold) '.txt'];
        fopen(path, 'w');
        writematrix(MI{index, 2, fold}, path, 'Delimiter', 'tab');
        fclose('all');
    end
end

% Storing the total mutual information
for fold = 1:number_folds
    path = ['results-MI-ICA/DataNewNoise/totalMI_it' int2str(fold) '.txt'];

    % Replacing the empty arrays in cell for zeros
    totalMIaux = totalMI(:,:,fold);
    emptycells = cellfun('isempty', totalMIaux);
    totalMIaux(emptycells) = {0};
    totalMIaux = cell2mat(totalMIaux);

    % Storing the data
    fopen(path, 'w');
    writematrix(totalMIaux, path, 'Delimiter', 'tab');
    fclose('all');
end



