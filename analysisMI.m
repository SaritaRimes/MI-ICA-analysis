% Calculo da informacao mutua %

clear all
close all
clc

addpath('FastICA_25', 'internalMethods');

MI = cell(11, 2); % column 1: no ICA, column 2: with ICA
totalMI = zeros(11, 2); % column 1: no ICA, column 2: with ICA

mPu = 50;
snr = 2;
nevents = 2000000;
occupancies = [10 30 50 80];

for oc = occupancies
    fprintf(['Mean of pile-up: ' int2str(mPu) ',\t Occupancy: ' int2str(oc) '\n']);

    data = load(['../../../RuidoSimuladoNovoSimulador/TileCal/ruido_media' int2str(mPu) '/ruido_ocup' int2str(oc) '_' ...
                 int2str(nevents) 'sinais.txt']); % load noise data
    data = data(1:100000,:);
    number_events = size(data,1);
    number_samples = size(data,2);
    number_intervals = 100; % number of discretization intervals
    h = (1-(-1))/number_intervals; % step of discretization
    index = (oc/10) + 1;

    % Applying ICA to the data
    [dataICA, A, W] = fastica(data'); % apply ICA
    dataICA = dataICA'; % transpoe a matriz para ter as variaveis nas colunas
    
    % Normalizing the variables between -1 and 1
    normalized_noise = normalizeSamples(data);  % normalized variables
    normalized_noise_ica = normalizeSamples(dataICA); % normalized variables
    
    % Discretizing the variables and getting the number of samples per interval
    [discretized_noise, number_samples_interval] ...
        = discretizeSamples(normalized_noise, number_intervals, oc, 0);
    [discretized_noise_ica, number_samples_interval_ica] ...
        = discretizeSamples(normalized_noise_ica, number_intervals, oc, 0);

    % Calculating the marginal probabilities
    marginal_probability = number_samples_interval/number_events;
    marginal_probability_ica = number_samples_interval_ica/number_events;
        
    % Checking whether the marginal probabilities sums 1
    checkProbabilities("marginal", marginal_probability);
    checkProbabilities("marginal", marginal_probability_ica);
        
    % Calculating the joint probabilities
    joint_probability = jointProbability(discretized_noise, ...
                                         number_intervals, oc, 0);
    joint_probability_ica = jointProbability(discretized_noise_ica, ...
                                             number_intervals, oc, 0);

    % Checking whether the joint probabilities sums 1
    checkProbabilities("joint", joint_probability);
    checkProbabilities("joint", joint_probability_ica);
    
    % Calculating the Mutual Information
    MI{index, 1} = mutualInformation(size(discretized_noise, 2), number_intervals, ...
                                     marginal_probability, joint_probability, oc, 0);
    MI{index, 2} = mutualInformation(size(discretized_noise_ica, 2), number_intervals, ...
                                     marginal_probability_ica, joint_probability_ica, oc, 0);
    
    % Calculating the crosstalk, which is the mutual information contained in all the process
    totalMI(index, 1) = mutualInformationCrosstalk(MI{index, 1}, number_samples);
    totalMI(index, 2) = mutualInformationCrosstalk(MI{index, 2}, number_samples);
end

% Plotting the graphs individually
% Before ICA
figure('Position', [100 100 800 500])
plot(occupancies, totalMI((occupancies/10 + 1), 1), 'Color', [0.6 0.6 0.6], 'Marker', '.', ...
     'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
xlim([0 100]);
title(['Total MI before ICA (mPu' int2str(mPu) ')']);
legend('MI before ICA', 'Location', 'northeast');
xlabel('Occupancy (%)');
ylabel('Total Mutual Information (ADC counts)');
% After ICA
figure('Position', [100 100 800 500])
plot(occupancies, totalMI((occupancies/10 + 1), 2), 'Color', [0.6 0.6 0.6], 'Marker', '.', ...
     'MarkerSize', 20, 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
xlim([0 100]);
title(['Total MI after ICA (mPu' int2str(mPu) ')']);
legend('MI after ICA', 'Location', 'northeast');
xlabel('Occupancy (%)');
ylabel('Total Mutual Information (ADC counts)');


% Plotting the graphs before and after the ICA
figure('Position', [100 100 800 500])
plot(occupancies, totalMI((occupancies/10 + 1), 1), 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
hold on
plot(occupancies, totalMI((occupancies/10 + 1), 2), 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
hold off
title(['Total MI before and after ICA (mPu' int2str(mPu) ')']);
legend({'before ICA', 'after ICA'}, 'Location', 'northeast');
xlabel('Occupancy (%)');
ylabel('Total Mutual Information (ADC counts)');
xlim([0, 100]);
ylim([0, 5*ceil(max(totalMI(:, 1))/5)]);


