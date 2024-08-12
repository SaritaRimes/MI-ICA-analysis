% Calculo da informacao mutua %

clear all
close all
clc

addpath('FastICA_25', 'internalMethods');

MI = cell(1,11);
totalMI = zeros(11,1);

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
    h = (1-(-1))/number_intervals;
    index = (oc/10) + 1;
    
    % Normalizing the variables between -1 and 1
    noICA = true;
    if ~noICA
        % Applying ICA to the data
        [dataICA, A, W] = fastica(data'); % apply ICA
        dataICA = dataICA'; % transpoe a matriz para ter as variaveis nas colunas

        normalized_noise = normalizeSamples(dataICA); % normalized variables
    else
        normalized_noise = normalizeSamples(data);  % normalized variables
    end
    
    % Discretizing the variables and getting the number of samples per interval
    [discretized_noise, number_samples_interval] = discretizeSamples(normalized_noise, number_intervals, oc);

    % Calculating the marginal probabilities
    marginal_probability = number_samples_interval/number_events;
        
    % Checking whether the marginal probabilities sums 1
    checkProbabilities("marginal", marginal_probability);
        
    % Calculating the joint probabilities
    joint_probability = jointProbability(discretized_noise, number_intervals, oc);

    % Checking whether the joint probabilities sums 1
    checkProbabilities("joint", joint_probability);
    
    % Calculating the Mutual Information
    MI{1, index} = mutualInformation(size(discretized_noise, 2), number_intervals, marginal_probability, ...
                                     joint_probability, oc);
    
    % Calculating the crosstalk, which is the mutual information contained in all the process
    totalMI(index) = mutualInformationCrosstalk(MI{1, index}, number_samples);
end

% Saving the data of MI in .txt files
% dlmwrite('results-mi-ica/totalMI_noICA.txt', totalMI);
% for oc = 10:10:100
%     ind = (oc/10) + 1;
%     dlmwrite(['results-MI-ICA/MI-occup' int2str(oc) '_noICA.txt'], MI{1, ind});
% end
% if ~noICA
%     dlmwrite('results-MI-ICA/totalMI_withICA.txt', totalMI);     
%     for oc = 10:10:100
%         ind = (oc/10) + 1;
%         dlmwrite(['resultados-mi-ica/MI-occup' int2str(oc) '_withICA.txt'], MI{1, ind});
%     end
% end

% Plotting the graphs individually
figure
plot(occupancies, totalMI((occupancies/10 + 1), 1), 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
xlim([0 100]);
title(['Total MI (mPu' int2str(mPu) ')']);
legend('MI before ICA', 'Location', 'northeast');
xlabel('Occupancy (%)');
ylabel('Total Mutual Information (ADC counts)');


% Plotting the graphs before and after the ICA
figure
plot(occupancies, MItotal_semICA((occupancies/10 + 1), 1), 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
hold on
plot(occupancies, MItotal_comICA((occupancies/10 + 1), 1), 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
hold off
title(['Total MI mPu' int2str(mPu) ' snr' int2str(snr)]);
legend({'before ICA', 'after ICA'}, 'Position', [0.78 0.7 0.1 0.2]);
xlabel('Occupancy (%)');
ylabel('Total Mutual Information (ADC counts)');
xlim([0, 110]);
ylim([min(MItotal_comICA) - 10, max(MItotal_semICA) + 10]);


