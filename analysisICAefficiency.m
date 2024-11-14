%% Effiiency Analysis %%
% Methods: OF, COF, Gaussian MLE e ICA MLE %
% Data: it was generated using the new simulator (github.com/ingoncalves/calorimetry-pulse-simulator) %

clear all
close all
% clc

addpath("FastICA_25");

mean_error_gauss = zeros(11, 1);
mean_error_of = zeros(11, 1);
mean_error_cof = zeros(11, 1);
mean_error_ica = zeros(11, 1);
std_error_gauss = zeros(11, 1);
std_error_of = zeros(11, 1);
std_error_cof = zeros(11, 1);
std_error_ica = zeros(11, 1);

mPu = 100;
snr = 3;
bins = 100;
number_dimensions = 3;

occupancies = [10 30 50 80];

for oc = occupancies

    noise = load(['../dadosRuido/comPedestal/mPu' int2str(mPu) '_snr' int2str(snr) ...
                  '/noise-ocup' int2str(oc) '.csv']); % load noise data
    pedestal = 50;
    %ruido = ruido - ped;
    
    % Dividing the before ICA data in two datasets
    div = cvpartition(size(noise,1), 'Holdout', 0.5); % choose 50% of signals randomly
    ind = div.test; % return the indexes of 50% of choosing signals
    noise_train = noise(ind,:);
    noise_test = noise(~ind,:);
    number_events = size(noise_test,1); %quantidade de sinais no conjunto de teste

    mean_noise_train = mean(noise_train);

    % Applying ICA to the noise data of training
    [noise_ica, A, W] = fastica(noise_train', 'numOfIC', number_dimensions); % ICA function
    noise_ica = noise_ica'; % variables must be in columns

    % Normalizing the histograms of noise after ICA
    hist_probabilities = -1*ones(bins, size(noise_ica,2));
    hist_bins = zeros(bins + 1, size(noise_ica,2));
    for i = 1:number_dimensions
        h = histogram(noise_ica(:,i), bins, 'Normalization', 'probability');
        hist_probabilities(:,i) = h.Values;
        hist_bins(:,i) = h.BinEdges;
        figure
    end

    % Finding the x coordinates of histograms
    hist_coordinate_x = zeros(size(hist_bins, 1) - 1, number_dimensions);
    spline_hist = zeros(number_dimensions);
    for j = 1:number_dimensions
        for i = 1:size(hist_bins, 1) - 1
            hist_coordinate_x(i,j) = (hist_bins(i,j) + hist_bins(i + 1,j))/2;
        end
        
        %Utilizando spline para interpolar
        %splineCoordenadaX = min(histCoordenadaX(:,j)):0.1:max(histCoordenadaX(:,j));
        spline_hist(j) = spline(hist_coordinate_x(:,j), hist_probabilities(:,j));
        %ppval(splineHist(3),min(histCoordenadaX(:,3)):0.1:max(histCoordenadaX(:,3)));
    end
    
    % Structures pre defined
    s = [0  0.0172  0.4524  1  0.5633  0.1493  0.0424]; %vetor de amostras do pulso de 
                                                        %referencia normalizado
    OF2 = [-0.3781  -0.3572  0.1808  0.8125  0.2767  -0.2056  -0.3292];
    
    % Applying the ICA to the normalized pulse
%     sICA = s*W;
    
    % Mounting the complete signal
    amplitude_true = exprnd(snr*mPu, number_events, 1);
    r = zeros(number_events, size(noise_test,2));
    for i = 1:number_events
        r(i,:) = amplitude_true(i)*pegaPulseJitter + noise_test(i,:); % complete signal in 7 dimensions
    end    
    
    % Estimating the amplitude using the linear methods
    covariance_gauss = cov(noise_train); % covariance matrix of training data
    OF = (inv(covariance_gauss)*s')/(s*inv(covariance_gauss)*s');
    amplitude_gauss = (r - pedestal)*OF;
    amplitude_of = r*OF2';
    amplitude_cof = aplicaCOF(r - pedestal,4.5);

    % Estimating the Gaussian PDF
    noise_gauss = r - amplitude_gauss*s;
    pdf_gauss = ones(number_events, 1);
    for i = 1:size(noise_gauss, 1)
        pdf_gauss(i) = (1/(sqrt(det(covariance_gauss))*(2*pi)^(3.5)))...
                          *exp(-.5*noise_gauss(i, :)*inv(covariance_gauss) ...
                               *noise_gauss(i, :)');
    end
    
    % Estimating the amplitude using MLE + ICA method
    amplitude_ica = amplitude_gauss;
    marginal_probability = zeros(1, number_dimensions);
    pdf_ica = ones(number_events, 1);
    for i = 1:number_events
        fprintf("Media do sinal = %d \nOcupacao = %d \nEvento = %d/%d \n", ...
                snr*mPu, oc, i, number_events);
        
        maximum_probability = -1;
        for amplitude_auxiliary = amplitude_gauss(i)-100:1:amplitude_gauss(i)+100
            
            ruidoAuxiliarTemporario = r(i,:) - amplitude_auxiliary*s;
            ruidoAuxiliar = (ruidoAuxiliarTemporario-mean_noise_train)*W';
            
            for j = 1:number_dimensions
                marginal_probability(j) = ppval(spline_hist(j), ruidoAuxiliar(j));
            end
            
            probability_ica = prod(marginal_probability);
            
            if (probability_ica > maximum_probability)
                amplitude_ica(i,1) = amplitude_auxiliary;
                maximum_probability = probability_ica;
            end
        end
        
        pdf_ica(i) = maximum_probability;
    end

    % Estimating the chi2 of each method
    chi2_gaussiano = zeros(number_events, 1);
    chi2_ica = zeros(number_events, 1);
    for i = 1:number_events
        % gaussian
        noise_temporary = r(i, :) - (amplitude_gauss(i)*s + pedestal);
        chi2_gaussiano(i) = sqrt(sum((noise_temporary.^2)./7));
        % ica
        noise_temporary = r(i, :) - (amplitude_ica(i)*s + pedestal);
        chi2_ica(i) = sqrt(sum((noise_temporary.^2)./7));
    end
    
    % Estimating the error of each method
    error_gauss(:,1) = amplitude_gauss - amplitude_true;
    error_of(:,1) = amplitude_of - amplitude_true;
    error_cof(:,1) = amplitude_cof - amplitude_true;
    error_ica(:,1) = amplitude_ica - amplitude_true;

    % Estimating the mean of the errors
    indice = oc/10 + 1;
    mean_error_gauss(indice, 1) = mean(error_gauss);
    mean_error_of(indice, 1) = mean(error_of);
    mean_error_cof(indice, 1) = mean(error_cof);
    mean_error_ica(indice, 1) = mean(error_ica);

    % Estimating the standard deviation of the errors
    std_error_gauss(indice, 1) = std(error_gauss);
    std_error_of(indice, 1) = std(error_of);
    std_error_cof(indice, 1) = std(error_cof);
    std_error_ica(indice, 1) = std(error_ica);
end

% Plotando os histogramas dos erros
histogram(error_gauss, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'LineWidth', 1.5);
hold on
histogram(error_of, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'g', 'LineWidth', 1.5);
hold on
histogram(error_cof, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on
histogram(error_ica, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'm', 'LineWidth', 1.5);
hold off
xlim([-500 600]);
legend({'MLE Gaussiano', 'OF', 'COF', 'MLE ICA'}, 'Position', [0.17 0.7 0.1 0.2]);
title(['Histograma dos erros com ' int2str(number_dimensions) ' dimensões']);

% Plotando os graficos de erro versus probabilidade
figure
scatter(error_gauss, pdf_gauss, 'MarkerEdgeColor', 'b', 'Marker', '.');
title(['MLE Gaussiano, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('Probabilidade');
figure
scatter(error_ica, pdf_ica, 'MarkerEdgeColor', 'r', 'Marker', '.');
title(['MLE + ICA, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('Probabilidade');

% Plotando os graficos de erro versus chi2
figure
scatter(error_gauss, chi2_gaussiano, 'MarkerEdgeColor', 'b', 'Marker', '.');
title(['MLE Gaussiano, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('\chi^2');
figure
scatter(error_ica, chi2_ica, 'MarkerEdgeColor', 'r', 'Marker', '.');
title(['MLE + ICA, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('\chi^2');

return;        

% Plotando os graficos do pulso normalizado antes e depois da ICA
x = 1:7;
plot(x, s, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
hold on
plot(x, sICA, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
hold off
%title(['Total MI mPu' int2str(mPu) ' snr' int2str(snr)]);
legend({'before ICA', 'after ICA'}, 'Position', [0.78 0.7 0.1 0.2]);
%xlabel('Occupancy (%)');
%ylabel('Total Mutual Information (ADC counts)');
xlim([0, 8]);

% Histograma das variaveis de ruido apos a aplicacao da ICA
figure
histogram(noise_ica(:,1), 100);
title(['Variavel 1, Ocupação ' int2str(oc)]);
figure
histogram(noise_ica(:,2), 100);
title(['Variavel 2, Ocupação ' int2str(oc)]);
figure
histogram(noise_ica(:,3), 100);
title(['Variavel 3, Ocupação ' int2str(oc)]);
figure
histogram(noise_ica(:,4), 100);
title(['Variavel 4, Ocupação ' int2str(oc)]);
figure
histogram(noise_ica(:,5), 100);
title(['Variavel 5, Ocupação ' int2str(oc)]);
figure
histogram(noise_ica(:,6), 100);
title(['Variavel 6, Ocupação ' int2str(oc)]);
figure
histogram(noise_ica(:,7), 100);
title(['Variavel 7, Ocupação ' int2str(oc)]);


% Plotando os graficos de media e desvio padrao
% x = 0:10:100;
% plot(x, mediaGaussiano, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
% hold on
% plot(x, mediaOF2, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
% hold on
% plot(x, mediaCOF, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.4940 0.1840 0.5560]);
% hold on
% plot(x, mediaLognormal, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.6350 0.0780 0.1840]);
% hold off
% title(['Média mPu' int2str(mPu) ' snr' int2str(snr)]);
% legend({'MLE Gaussiano', 'OF', 'COF', 'MLE Lognormal'}, 'Position', [0.17 0.7 0.1 0.2]);
% xlabel('Ocupação (%)');
% ylabel('Média do erro (ADC counts)');

% plot(x,desvioPadraoGaussiano, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
% hold on
% plot(x, desvioPadraoOF2, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
% hold on
% plot(x, desvioPadraoCOF, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.4940 0.1840 0.5560]);
% hold on
% plot(x, desvioPadraoLognormal, 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0.6350 0.0780 0.1840]);
% hold off
% title(['Desvio padrão mPu' int2str(mPu) ' snr' int2str(snr)]);
% legend({'MLE Gaussiano', 'OF', 'COF', 'MLE Lognormal'}, 'Position', [0.17 0.7 0.1 0.2]);
% xlabel('Ocupação (%)');
% ylabel('Desvio padrão do erro (ADC counts)');
