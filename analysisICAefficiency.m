%% Analise de eficiencia %%
% Metodos: OF, COF, MLE Gaussiano e MLE Lognormal %
% Dados: gerados com o simulador em Python %

clear all
close all
% clc

addpath("FastICA_25");

%Iniciando algumas estruturas que serao utilizadas
mediaGaussiano = zeros(11, 1);
mediaOF2 = zeros(11, 1);
mediaCOF = zeros(11, 1);
mediaICA = zeros(11, 1);
desvioPadraoGaussiano = zeros(11, 1);
desvioPadraoOF2 = zeros(11, 1);
desvioPadraoCOF = zeros(11, 1);
desvioPadraoICA = zeros(11, 1);

%Definindo algumas variaveis
mPu = 100;
snr = 3;
bins = 100;
numeroDimensoes = 3;

for oc = 20:10:20

    % Carregando os dados de ruido e definindo o pedestal
    ruido = load(['../dadosRuido/comPedestal/mPu' int2str(mPu) '_snr' int2str(snr) '/noise-ocup' int2str(oc) '.csv']); %carrega os dados de ruido
    ped = 50;
    %ruido = ruido - ped;
    
    % Fazendo a divisao dos dados, antes da ICA, em dois conjuntos
    div = cvpartition(size(ruido,1), 'Holdout', 0.5); %escolhe 50% dos sinais aleatoriamente
    ind = div.test; %retorna os indices dos 50% escolhidos
    ruidoTreino = ruido(ind,:); %os 50% sao selecionados para o conjunto de treino
    ruidoTeste = ruido(~ind,:); %os outros 50% sao selecionados para o conjunto de teste
    numeroEventos = size(ruidoTeste,1); %quantidade de sinais no conjunto de teste

    mediaRuidoTreino = mean(ruidoTreino);
    % Aplicando a ICA aos dados de ruido de treinamento
    [ruidoICA, A, W] = fastica(ruidoTreino', 'numOfIC', numeroDimensoes); %aplica a ICA
    ruidoICA = ruidoICA'; %transpoe a matriz para ter as variaveis nas colunas
%     return;
    % Normalizando os histogramas do ruido apos a ICA
    histProbabilidades = -1*ones(bins,size(ruidoICA,2));
    histIntervalosBins = zeros(bins + 1,size(ruidoICA,2));
    for i = 1:numeroDimensoes
        h = histogram(ruidoICA(:,i),bins,'Normalization','probability');
        histProbabilidades(:,i) = h.Values;
        histIntervalosBins(:,i) = h.BinEdges;
%         return;
        figure;
    end
%     return;
    % Encontrando coordenadas x dos histogramas
    histCoordenadaX = zeros(size(histIntervalosBins, 1) - 1, numeroDimensoes);
    for j = 1:numeroDimensoes
        for i = 1:size(histIntervalosBins, 1) - 1
            histCoordenadaX(i,j) = (histIntervalosBins(i,j) + histIntervalosBins(i + 1,j))/2;
        end
        
        %Utilizando spline para interpolar
        %splineCoordenadaX = min(histCoordenadaX(:,j)):0.1:max(histCoordenadaX(:,j));
        splineHist(j) = spline(histCoordenadaX(:,j), histProbabilidades(:,j));
        %ppval(splineHist(3),min(histCoordenadaX(:,3)):0.1:max(histCoordenadaX(:,3)));
    end
    
%     return;
    
    % Vetores de entrada
    s = [0  0.0172  0.4524  1  0.5633  0.1493  0.0424]; %vetor de amostras do pulso de 
                                                        %referencia normalizado
    OF2 = [-0.3781  -0.3572  0.1808  0.8125  0.2767  -0.2056  -0.3292];
    
    % Aplicando a ICA ao pulso normalizado
    %sICA = s*W;
    
    % Montando o sinal completo
    amplitudeVerdadeira = exprnd(snr*mPu, numeroEventos, 1);
    r = zeros(numeroEventos, size(ruidoTeste,2));
    for i = 1:numeroEventos
        r(i,:) = amplitudeVerdadeira(i)*pegaPulseJitter + ruidoTeste(i,:); %sinal completo em 7 dimensoes
    end    
    
    % Calculando a amplitude pelos metodos lineares
    covarianciaGaussiano = cov(ruidoTreino); %matriz de covariancia do ruido
    OF = (inv(covarianciaGaussiano)*s')/(s*inv(covarianciaGaussiano)*s');
    amplitudeGaussiano = (r - ped)*OF;
    amplitudeOF2 = r*OF2';
    amplitudeCOF = aplicaCOF(r - ped,4.5);

    % Calculando a PDF gaussiana
    ruidoGaussiano = r - amplitudeGaussiano*s;
    pdfGaussiano = ones(numeroEventos, 1);
    for i = 1:size(ruidoGaussiano, 1)
        pdfGaussiano(i) = (1/(sqrt(det(covarianciaGaussiano))*(2*pi)^(3.5)))...
                          *exp(-.5*ruidoGaussiano(i, :)*inv(covarianciaGaussiano) ...
                               *ruidoGaussiano(i, :)');
    end
    
    % Calculando a amplitude pelo metodo MLE com ICA
    amplitudeICA = amplitudeGaussiano;
    probMarginal = zeros(1, numeroDimensoes);
    pdfICA = ones(numeroEventos, 1);
    for i = 1:numeroEventos
        %fprintf("Aplicando o metodo MLE Lognormal.\n");
        fprintf("Media do sinal = %d \nOcupacao = %d \nEvento = %d/%d \n", snr*mPu, oc, i, numeroEventos);
        
        probMaxima = -1;
        for amplitudeAuxiliar = amplitudeGaussiano(i)-100:1:amplitudeGaussiano(i)+100
            
            ruidoAuxiliarTemporario = r(i,:) - amplitudeAuxiliar*s;
            ruidoAuxiliar = (ruidoAuxiliarTemporario-mediaRuidoTreino)*W';
            
            for j = 1:numeroDimensoes
                probMarginal(j) = ppval(splineHist(j), ruidoAuxiliar(j));
            end
            
%             return;
            
%             for j = 1:(size(histIntervalosBins, 1) - 1)
%                 for k = 1:size(ruidoAuxiliar, 2)
%                     if ruidoAuxiliar(k) >= histIntervalosBins(j,k) && ruidoAuxiliar(k) < histIntervalosBins(j + 1,k)
%                         probMarginal(k) = histProbabilidades(j,k);
%                     end
%                 end
%             end
            
            probICA = prod(probMarginal);
            
            if (probICA > probMaxima)
                amplitudeICA(i,1) = amplitudeAuxiliar;
                probMaxima = probICA;
            end
        end
        
        pdfICA(i) = probMaxima;
    end

    % Calculando o chi2 de cada metodo
    chi2Gaussiano = zeros(numeroEventos, 1);
    chi2ICA = zeros(numeroEventos, 1);
    for i = 1:numeroEventos
        ruidoTemporario = r(i, :) - (amplitudeGaussiano(i)*s + ped);
        chi2Gaussiano(i) = sqrt(sum((ruidoTemporario.^2)./7));
        ruidoTemporario = r(i, :) - (amplitudeICA(i)*s + ped);
        chi2ICA(i) = sqrt(sum((ruidoTemporario.^2)./7));
    end
    
    % Calculando os erros de cada metodo
    erroGaussiano(:,1) = amplitudeGaussiano - amplitudeVerdadeira;
    erroOF2(:,1) = amplitudeOF2 - amplitudeVerdadeira;
    erroCOF(:,1) = amplitudeCOF - amplitudeVerdadeira;
    erroICA(:,1) = amplitudeICA - amplitudeVerdadeira;

    % Calculando a media dos erros
    indice = oc/10 + 1;
    mediaGaussiano(indice, 1) = mean(erroGaussiano);
    mediaOF2(indice, 1) = mean(erroOF2);
    mediaCOF(indice, 1) = mean(erroCOF);
    mediaICA(indice, 1) = mean(erroICA);

    % Calculando o desvio padrao dos erros
    desvioPadraoGaussiano(indice, 1) = std(erroGaussiano);
    desvioPadraoOF2(indice, 1) = std(erroOF2);
    desvioPadraoCOF(indice, 1) = std(erroCOF);
    desvioPadraoICA(indice, 1) = std(erroICA);

end

% Plotando os histogramas dos erros
histogram(erroGaussiano, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'LineWidth', 1.5);
hold on
histogram(erroOF2, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'g', 'LineWidth', 1.5);
hold on
histogram(erroCOF, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on
histogram(erroICA, 100, 'DisplayStyle', 'stairs', 'EdgeColor', 'm', 'LineWidth', 1.5);
hold off
xlim([-500 600]);
legend({'MLE Gaussiano', 'OF', 'COF', 'MLE ICA'}, 'Position', [0.17 0.7 0.1 0.2]);
title(['Histograma dos erros com ' int2str(numeroDimensoes) ' dimensões']);

% Plotando os graficos de erro versus probabilidade
figure
scatter(erroGaussiano, pdfGaussiano, 'MarkerEdgeColor', 'b', 'Marker', '.');
title(['MLE Gaussiano, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('Probabilidade');
figure
scatter(erroICA, pdfICA, 'MarkerEdgeColor', 'r', 'Marker', '.');
title(['MLE + ICA, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('Probabilidade');

% Plotando os graficos de erro versus chi2
figure
scatter(erroGaussiano, chi2Gaussiano, 'MarkerEdgeColor', 'b', 'Marker', '.');
title(['MLE Gaussiano, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('\chi^2');
figure
scatter(erroICA, chi2ICA, 'MarkerEdgeColor', 'r', 'Marker', '.');
title(['MLE + ICA, ocupação ' int2str(oc) '%'], 'FontSize', 13);
xlabel('Erro (contagens de ADC)');
ylabel('\chi^2');

return;

    % Montando o sinal completo apos a ICA
%     amplitudeVerdadeiraICA = exprnd(mediaSinalICA, numeroEventosICA, 1);
%     rICA = zeros(numeroEventosICA, size(ruidoTesteICA,2));
%     for i = 1:numeroEventosICA
%         rICA(i,:) = amplitudeVerdadeiraICA(i)*pegaPulseJitter + ruidoTesteICA(i,:); %sinal completo em 7 dimensoes
%     end
%     rICA = r*W;
    
    % Calculando os parametros da distribuicao Lognormal
%     k = 1; %contador
%     mediaLogn = zeros(1,size(ruidoICA,2));
%     desvioPadraoLogn = zeros(1,size(ruidoICA,2));
%     for i = 1:size(ruidoICA,2)
%         k = ruidoICA(:,i) > 0;
%         mediaLogn(i) = mean(log(ruidoICA(k,i)));
%         desvioPadraoLogn(i) = std(log(ruidoICA(k,i)));
%     end

    % Calculando a media do sinal de interesse apos a ICA
%     tamanhoRuidoPositivoICA = -1;
%     for i = 1:size(ruidoICA, 2)
%         k = ruidoICA(:,i) > 0; %seleciona apenas as amostras > 0 da coluna i
%         ruidoPositivoTemporarioICA = ruidoICA(k,i);        
%         
%         %Verificando qual coluna sera usada para calcular a media
%         if size(ruidoPositivoTemporarioICA,1) > tamanhoRuidoPositivoICA
%             tamanhoRuidoPositivoICA = size(ruidoPositivoTemporarioICA,1);
%             ruidoPositivoICA = ruidoPositivoTemporarioICA;
%             maiorDimensao = i;
%         end
%     end
%     mediaSinalICA = mean(ruidoPositivoICA);

    % Calculando a media do sinal de interesse antes da ICA
%     tamanhoRuidoPositivo = -1;
%     for i = 1:size(ruido, 2)
%         k = ruido(:,i) > 0; %seleciona apenas as amostras > 0 da coluna i
%         ruidoPositivoTemporario = ruido(k,i);
%         
%         %Verificando qual coluna sera usada para calcular a media
%         if size(ruidoPositivoTemporario,1) > tamanhoRuidoPositivo
%             tamanhoRuidoPositivo = size(ruidoPositivoTemporario,1);
%             ruidoPositivo = ruidoPositivoTemporario;
%             maiorDimensao = i;
%         end
%     end
%     mediaSinal = mean(ruidoPositivo);

    % Fazendo a divisao dos dados, apos a ICA, em dois conjuntos
%     div = cvpartition(size(ruidoICA,1), 'Holdout', 0.5); %escolhe 50% dos sinais aleatoriamente
%     ind = div.test; %retorna os indices dos 50% escolhidos
%     ruidoTreinoICA = ruidoICA(ind,:); %os 50% sao selecionados para o conjunto de treino
%     ruidoTesteICA = ruidoICA(~ind,:); %os outros 50% sao selecionados para o conjunto de teste
%     numeroEventosICA = size(ruidoTesteICA,1); %quantidade de sinais no conjunto de teste
        

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
hist(ruidoICA(:,1), 100);
title(['Variavel 1, Ocupação ' int2str(oc)]);
figure
hist(ruidoICA(:,2), 100);
title(['Variavel 2, Ocupação ' int2str(oc)]);
figure
hist(ruidoICA(:,3), 100);
title(['Variavel 3, Ocupação ' int2str(oc)]);
figure
hist(ruidoICA(:,4), 100);
title(['Variavel 4, Ocupação ' int2str(oc)]);
figure
hist(ruidoICA(:,5), 100);
title(['Variavel 5, Ocupação ' int2str(oc)]);
figure
hist(ruidoICA(:,6), 100);
title(['Variavel 6, Ocupação ' int2str(oc)]);
figure
hist(ruidoICA(:,7), 100);
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











