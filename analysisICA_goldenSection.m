%implementacao do MLE lognormal força bruta

clear
close all
clc

% Setting some constants
mean_pileup = 50;
snr = 2;
number_bins = 100;
bins_width = 5;
occupancies = [10 30 50 80];
quantity_positive_samples = zeros(4, 1);
number_iterations = 10;
quantity_signals = 2000000;
number_dimensions_ica = 7;

% Setting some structures
errors = cell(5, 4, number_iterations); % in the rows the methods and in the columns the occupations
amplitudes = cell(7, 5, number_iterations); % in the rows the methods and in the columns the occupations
kl_divergence = zeros(5, 4); % in the rows the methods and in the columns the occupations
ks_statistic = zeros(5, 4); % in the rows the methods and in the columns the occupations

kl_divergence_gauss = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
kl_divergence_of = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
kl_divergence_cof = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
kl_divergence_logn = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
ks_statistic_gauss = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
ks_statistic_of = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
ks_statistic_cof = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations
ks_statistic_logn = zeros(5, number_iterations); % in the rows the occupations and in the columns the iterations

probabilities = cell(6, 4, number_iterations); % in the rows the methods and in the columns the occupations
chi2 = cell(5, 4, number_iterations); % in the rows the methods and in the columns the occupations
std_error_cuts_prob = cell(5, 4, number_iterations); % in the rows the methods and in the columns the occupations
std_error_cuts_chi2 = cell(5, 4, number_iterations); % in the rows the methods and in the columns the occupations


addpath('internalMethods', 'externalMethods', 'estimationMethods', 'dataQuality', 'FastICA_25');

% Calorimeter reference pulse
s = [0 .0172 .4524 1 .5633 .1493 .0424];
% OF2
OF2 = [ -0.3781   -0.3572    0.1808    0.8125    0.2767   -0.2056   -0.3292 ];

index_occupancy = 1;
for oc = occupancies
     
    % Loading the noise and setting the pedestal
    noise_total = load(['D:/Documentos/UERJ/Doutorado/Simulacoes/RuidoSimuladoNovoSimulador/TileCal/ruido_media' ...
                       int2str(mean_pileup) '/ruido_ocup' int2str(oc) '_' int2str(quantity_signals) 'sinais.txt']);
    ped = 50;

    % Calculating the number of signals in each set
    number_signals = quantity_signals / number_iterations;
    
    for it = 1:number_iterations

        % Selecting a part of the noise set
        start = (it - 1) * number_signals + 1;
        final = it * number_signals;
        noise = noise_total(start:final, :);
        
        % Dividing the noise set into training and testing
        div = cvpartition(size(noise, 1), 'HoldOut', 0.50);
        ind = div.test;
        noise_training = noise(ind,:);
        noise_test = noise(~ind,:);
        number_events = size(noise_test, 1);
        
        % Removing signals containing negative samples
        indexes_noise_positive = -1*ones(size(noise_training, 1), 1);
        for i = 1:size(noise_training, 1)
            if ~any(noise_training(i,:) <= 0)
                indexes_noise_positive(i) = i;
            end
        end
        indexes_noise_positive = indexes_noise_positive(indexes_noise_positive > 0);
        noise_training_positive = noise_training(indexes_noise_positive, :);
        
        % MLE parameters
        mean_logn = mean(log(noise_training_positive));
        mean_gauss = mean(noise_training);
        covariance_logn = cov(log(noise_training_positive));
        covariance_gauss = cov(noise_training);
        
        % Assembling the complete signal
        amplitude_true = exprnd(snr*mean_pileup, number_events, 1);
        r = ones(size(noise_test, 1), size(noise_test, 2));
        for i = 1:number_events
            r(i,:) = amplitude_true(i)*s + noise_test(i,:);
        end
        
        % Linear methods
        OF = (covariance_gauss\s')/((s/covariance_gauss)*s');
        OF2 = tile_of2(noise_training, 1);
        amplitude_gauss = (r - ped)*OF;
        amplitude_of = r*OF2';
        amplitude_cof = aplicaCOF(r - ped, 4.5);
        amplitudes_signals_cof = aplicaCOFAll(r - ped, 4.5);
        
        % Lognormal MLE
        amplitude_logn = ones(size(r, 1)  , 1);
        for i = 1:number_events 
          
            fprintf(['MLE Lognormal: ' ...
                     'Processing event ' int2str(i) '/' int2str(number_events) '...' ...
                     '\nOccupancy ' int2str(oc), ', iteration ' int2str(it) '\n']);

            % Verificando se ha amostras negativas no sinal recebido
            if any(r(i, :) <= 0)
                amplitude_logn(i) = amplitude_gauss(i);
                continue;
            end
    
            % Estimando a amplitude via razão áurea
            amplitude_logn(i) ...
                = golden_section(@(A)pdfLognormal(mean_logn, ...
                                                  covariance_logn, ...
                                                  r(i,:), s, A), ...
                                 0, 1023, 1);

            if amplitude_logn(i) == 0
                amplitude_logn(i) = amplitude_gauss(i);
            end

        end

        % ICA
        % Applying ICA to the noise data of training
        [noise_ica, A, W] = fastica(noise_training', 'numOfIC', number_dimensions_ica); % ICA function
        noise_ica = noise_ica'; % variables must be in columns
        % Estimating the amplitude and the probability
        [amplitude_ica, pdf_ica] = IcaEstimator(r, s, noise_ica, mean_gauss, W, number_dimensions_ica, ...
                                                number_bins, amplitude_gauss, oc, it);        
        
        % Calculating the errors
        error_gauss = amplitude_gauss - amplitude_true;
        error_of = amplitude_of - amplitude_true;
        error_cof = amplitude_cof - amplitude_true;
        error_logn = amplitude_logn - amplitude_true;
        error_ica = amplitude_ica - amplitude_true;
        
        % Storing the errors
        errors{1, index_occupancy, it} = error_gauss;
        errors{2, index_occupancy, it} = error_of;
        errors{3, index_occupancy, it} = error_cof;
        errors{4, index_occupancy, it} = error_logn;
        errors{5, index_occupancy, it} = error_ica;

        % Calculating probabilities, chi2 and the standard deviation of data quality
        [probabilities(:, index_occupancy, it), ...
         chi2(:, index_occupancy, it), ...
         std_error_cuts_prob(:, index_occupancy, it), ...
         std_error_cuts_chi2(:, index_occupancy, it)] ...
            = analysisDataQuality({amplitude_of, amplitude_gauss, ...
                                   amplitudes_signals_cof, amplitude_logn, ...
                                   amplitude_true}, r, ...
                                  mean_logn, mean_gauss, ...
                                  covariance_logn, covariance_gauss, ...
                                  ped, oc);

        % Storing the amplitudes
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'efficiency/gaussiano/ocup' int2str(oc) '/amplitude_gaussiano_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', amplitude_gauss);
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'efficiency/of/ocup' int2str(oc) '/amplitude_of_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', amplitude_of);
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'efficiency/cof/ocup' int2str(oc) '/amplitude_cof_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', amplitude_cof);
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'efficiency/lognormal/ocup' int2str(oc) '/amplitude_lognormal_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', amplitude_logn);
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'efficiency/verdadeira/ocup' int2str(oc) '/amplitude_verdadeira_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.d \t %.13f\n', [find(~ind == true)'; amplitude_true']);
%         amplitudes{1, index_occupancy, it} = amplitude_gauss;
%         amplitudes{2, index_occupancy, it} = amplitude_of;
%         amplitudes{3, index_occupancy, it} = amplitude_cof;
%         amplitudes{4, index_occupancy, it} = amplitude_logn;
%         amplitudes{5, index_occupancy, it} = amplitude_true;

        % Storing the errors for data quality estimation
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/erro_data_quality/erro_mle_gaussiano_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', errors{1, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/erro_data_quality/erro_of_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', errors{2, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/erro_data_quality/erro_cof_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', errors{3, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/erro_data_quality/erro_mle_lognormal_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', errors{4, index_occupancy, it});
% 
%         % Storing the probabilities
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/probabilidade/prob_mle_gaussiano_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.35f\n', probabilities{1, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/probabilidade/prob_of_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.35f\n', probabilities{2, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/probabilidade/prob_cof_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.35f\n', probabilities{3, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/probabilidade/prob_mle_lognormal_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.35f\n', probabilities{4, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/probabilidade/prob_mle_lognormalgauss_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.35f\n', probabilities{5, index_occupancy, it});

        % Storing the chi2
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/chi2/chi2_mle_gaussiano_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', chi2{1, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/chi2/chi2_of_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', chi2{2, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/chi2/chi2_cof_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', chi2{3, index_occupancy, it});
%         path = ['D:/Documentos/UERJ/Doutorado/ArtigoIEEE/dados_txt/cross_validation/' ...
%                 'data_quality/ocup' int2str(oc) '/chi2/chi2_mle_lognormal_ocup' ...
%                 int2str(oc) '_it' int2str(it) '.txt'];
%         fprintf(fopen(path, 'w'), '%.13f\n', chi2{4, index_occupancy, it});

        % Plotting the histograms of errors
%         figure('Position', [50 100 900 600])
%         number_bins = round((max(erro_gauss) - min(erro_gauss)) / bins_width);
%         histogram(error_gauss, number_bins, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'LineWidth', 1.5);
%         hold on
%         number_bins = round((max(erro_of) - min(erro_of)) / bins_width);
%         histogram(error_of, number_bins, 'DisplayStyle', 'stairs', 'EdgeColor', 'g', 'LineWidth', 1.5);
%         hold on
%         number_bins = round((max(erro_cof) - min(erro_cof)) / bins_width);
%         histogram(error_cof, number_bins, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
%         hold on
%         number_bins = round((max(erro_logn) - min(erro_logn)) / bins_width);
%         histogram(error_logn, number_bins, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'LineWidth', 1.5);
%         hold on
%         histogram(error_ica, number_bins, 'DisplayStyle', 'stairs', 'EdgeColor', 'm', 'LineWidth', 1.5);
%         hold off
%         xlim([-300 600]);
%         legend({'MLE Gaussiano', 'OF', 'COF', 'MLE Lognormal', 'MLE + ICA'}, ...
%                'FontSize', 12, 'Location', 'northeast');
%         title(['Histogram of error with occupancy ' int2str(oc) '%'], 'FontSize', 13);
%         xlabel('Error (contagens de ADC)');
%         ylabel('Number of events');

    end

    index_occupancy = index_occupancy + 1;
end

return

% Plotting the quantity of positive samples
figure('Position', [50 100 900 600]);
plot(occupancies, quantity_positive_samples, '.-', 'MarkerSize', 20);
xlim([0 100]);
ylim([0 120]);
xlabel('Occupancy (%)');
ylabel('Positive samples (%)');

% Plotting the kl divergence for the occupancies
figure('Position', [50 100 900 600]);
plot(occupancies, kl_divergence(1, :), '.--', 'MarkerSize', 20, 'Color', 'b');
hold on
plot(occupancies, kl_divergence(2, :), '.--', 'MarkerSize', 20, 'Color', 'g');
hold on
plot(occupancies, kl_divergence(3, :), '.--', 'MarkerSize', 20, 'Color', 'r');
hold on
plot(occupancies, kl_divergence(4, :), '.--', 'MarkerSize', 20, 'Color', 'm');
legend({'MLE Gaussiano', 'OF', 'COF', 'MLE Lognormal'}, 'FontSize', 15, 'Location', 'northwest');
xlabel('Occupancy (%)', 'FontSize', 15);
ylabel('KL Divergence', 'FontSize', 15);

% Plotting the ks statistics for the occupancies
figure('Position', [50 100 900 600]);
plot(occupancies, ks_statistic(1, :), '.--', 'MarkerSize', 20, 'Color', 'b');
hold on
plot(occupancies, ks_statistic(2, :), '.--', 'MarkerSize', 20, 'Color', 'g');
hold on
plot(occupancies, ks_statistic(3, :), '.--', 'MarkerSize', 20, 'Color', 'r');
hold on
plot(occupancies, ks_statistic(4, :), '.--', 'MarkerSize', 20, 'Color', 'm');
legend({'MLE Gaussiano', 'OF', 'COF', 'MLE Lognormal'}, 'FontSize', 15, 'Location', 'northwest');
xlabel('Occupancy (%)', 'FontSize', 15);
ylabel('KS Statistic', 'FontSize', 15);

% Plotando a divergencia KL nas iteracoes
%%% Gaussian MLE
figure('Position', [50 100 900 600]);
plot(xIteracoes, kl_divergence_gauss(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_gauss(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_gauss(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, kl_divergence_gauss(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(kl_divergence_gauss)) - 0.01) (max(max(kl_divergence_gauss)) + 0.01)]);
title('Gaussian MLE', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KL Divergence', 'FontSize', 15);
%%% OF2
figure('Position', [50 100 900 600]);
plot(xIteracoes, kl_divergence_of(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_of(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_of(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, kl_divergence_of(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(kl_divergence_of)) - 0.01) (max(max(kl_divergence_of)) + 0.01)]);
title('OF2', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KL Divergence', 'FontSize', 15);
%%% COF
figure('Position', [50 100 900 600]);
plot(xIteracoes, kl_divergence_cof(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_cof(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_cof(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, kl_divergence_cof(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(kl_divergence_cof)) - 0.01) (max(max(kl_divergence_cof)) + 0.01)]);
title('COF', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KL Divergence', 'FontSize', 15);
%%% Lognormal MLE
figure('Position', [50 100 900 600]);
plot(xIteracoes, kl_divergence_logn(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_logn(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, kl_divergence_logn(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, kl_divergence_logn(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(kl_divergence_logn)) - 0.005) (max(max(kl_divergence_logn)) + 0.005)]);
title('Lognormal MLE', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KL Divergence', 'FontSize', 15);

% Plotting the ks statistic for the iterations
xIteracoes = 1:number_iterations;
%%% Gaussian MLE
figure('Position', [50 100 900 600]);
plot(xIteracoes, ks_statistic_gauss(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_gauss(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_gauss(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, ks_statistic_gauss(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(ks_statistic_gauss)) - 0.01) (max(max(ks_statistic_gauss)) + 0.01)]);
title('Gaussian MLE', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KS Statistic', 'FontSize', 15);
%%% OF2
figure('Position', [50 100 900 600]);
plot(xIteracoes, ks_statistic_of(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_of(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_of(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, ks_statistic_of(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(ks_statistic_of)) - 0.01) (max(max(ks_statistic_of)) + 0.01)]);
title('OF2', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KS Statistic', 'FontSize', 15);
%%% COF
figure('Position', [50 100 900 600]);
plot(xIteracoes, ks_statistic_cof(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_cof(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_cof(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, ks_statistic_cof(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(ks_statistic_cof)) - 0.01) (max(max(ks_statistic_cof)) + 0.01)]);
title('COF', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KS Statistic', 'FontSize', 15);
%%% Lognormal MLE
figure('Position', [50 100 900 600]);
plot(xIteracoes, ks_statistic_logn(1, :), '.--', 'Color', 'b', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_logn(2, :), '.--', 'Color', 'r', 'MarkerSize', 20);
hold on
plot(xIteracoes, ks_statistic_logn(3, :), '.--', 'Color', 'g', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
hold on
plot(xIteracoes, ks_statistic_logn(4, :), '.--', 'Color', 'm', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
legend({'10%', '30%', '50%', '80%'}, 'FontSize', 15, 'Location', 'best');
xlim([0 11]);
ylim([(min(min(ks_statistic_logn)) - 0.005) (max(max(ks_statistic_logn)) + 0.005)]);
title('Lognormal MLE', 'FontSize', 15);
xlabel('Iterations', 'FontSize', 15);
ylabel('KS Statistic', 'FontSize', 15);

