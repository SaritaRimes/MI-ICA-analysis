% Calculo da informacao mutua %

clear all
close all
clc

addpath('FastICA_25', 'internalMethods');

MI = cell(1,11);
MItotal = zeros(11,1);

mPu = 50;
snr = 2;
nevents = 2000000;
occupancies = [10 30 50 80];

for oc = occupancies
    fprintf(['Mean of pile-up: ' int2str(mPu) ',\t Occupancy: ' int2str(oc) '\n']);

    dados = load(['../../../RuidoSimuladoNovoSimulador/TileCal/ruido_media' int2str(mPu) '/ruido_ocup' int2str(oc) '_' ...
                  int2str(nevents) 'sinais.txt']); % carrega os dados de ruido
    dados = dados(1:100000,:);
    number_events = size(dados,1);
    number_samples = size(dados,2);
    number_intervals = 100; % numero de intervalos de discretizacao
    h = (1-(-1))/number_intervals;
    index = (oc/10) + 1;
    
    %Aplicando a ICA aos dados
%     [dadosICA, A, W] = fastica(dados'); % aplica a ICA
%     dadosICA = dadosICA'; % transpoe a matriz para ter as variaveis nas colunas
    
    %Normalizando as variaveis entre -1 e 1 com ICA
%     semICA = false;
%     ruido = zeros(neventos,tam);
%     for i = 1:tam
%         nume = 2*(dadosICA(:,i) - min(dadosICA(:,i)));
%         deno = max(dadosICA(:,i)) - min(dadosICA(:,i));
%         ruido(:,i) = (nume/deno) - 1;
%     end
%     
%     %Normalizando as variaveis entre -1 e 1 sem ICA
%     semICA = true;
%     ruido = zeros(neventos,tam);
%     for i = 1:tam
%         nume = 2*(dados(:,i) - min(dados(:,i)));
%         deno = max(dados(:,i)) - min(dados(:,i));
%         ruido(:,i) = (nume/deno) - 1;
%     end

    semICA = true;
    ruido = zeros(number_events,number_samples);
    if ~semICA
        % Aplicando a ICA aos dados
        [dadosICA, A, W] = fastica(dados'); % aplica a ICA
        dadosICA = dadosICA'; % transpoe a matriz para ter as variaveis nas colunas

        % Normalizando as variaveis entre -1 e 1
        for i = 1:number_samples
            nume = 2*(dadosICA(:,i) - min(dadosICA(:,i)));
            deno = max(dadosICA(:,i)) - min(dadosICA(:,i));
            ruido(:,i) = (nume/deno) - 1;
        end
    else
        % Normalizando as variaveis entre -1 e 1
        for i = 1:number_samples
            nume = 2*(dados(:,i) - min(dados(:,i)));
            deno = max(dados(:,i)) - min(dados(:,i));
            ruido(:,i) = (nume/deno) - 1;
        end        
    end
    
    % Salvando as variaveis discretizadas e calculando as prob marginais
    [discretized_noise, number_samples_interval] = discretizeSamples(ruido, number_intervals, oc);
    marginal_probability = number_samples_interval/number_events;
%     ruido_discretizado = -1*ones(size(ruido, 1), size(ruido, 2)); %serao as amostras discretizadas
%     classe = zeros(number_intervals,tam);
%     for amostra = 1:tam
%         inte = 1; % contabiliza os intervalos
% 
%         for i = -1:h:(1-h)
%             fprintf(['# Discretization \nOccupancy: ' int2str(oc) ',\t Sample: ' int2str(amostra)  ',\t Interval: ' num2str(i) '\n']);
%             if i == 1-h
%                 k = find((ruido(:,amostra) >= i) & (ruido(:,amostra) <= (i + h)));
%                 ruido_discretizado(k, amostra) = inte;
%             else
%                 k = find((ruido(:,amostra) >= i) & (ruido(:,amostra) < (i + h)));
%                 ruido_discretizado(k, amostra) = inte;
%             end
%             classe(inte, amostra) = size(k,1);
%             inte = inte + 1;
%             clear k
%         end
%     end
%     probmarg = classe/number_events;
        
    % Conferindo se as probabilidades marginais somam 1
    sum_marginal_probability = sum(marginal_probability);
    if round(sum_marginal_probability, 4) ~= 1
        throw('Sum of marginal probabilities is not equal to 1.');
    end
        
    % Calculando as probabilidades conjuntas
    joint_probability = jointProbability(discretized_noise, number_intervals, oc);
%     probconj = cell(tam,tam);
%     for am1 = 1:tam
%         for am2 = 1:tam
%             probconjaux = zeros(number_intervals,number_intervals);
%             for x = 1:number_intervals
%                 for y = 1:number_intervals
%                     for i = 1:number_events
% %                         fprintf('probabilidade conjunta \nocupação %d \namostra1 %d \n', oc, am1);
% %                         fprintf('amostra2 %d \nx = %d \ny = %d \nevento %d \n', am2, x, y, i);
%                         if (ruido_discretizado(i,am1) == x && ruido_discretizado(i,am2) == y)
%                             probconjaux(x,y) = probconjaux(x,y) + 1;
%                         end
%                     end
%                 end
%             end
%             probconj{am1,am2} = probconjaux/number_events;
%         end
%     end
    
    % Conferindo se as probabilidades conjuntas somam 1
    sum_joint_probability = zeros(number_samples, number_samples);
    for am1 = 1:number_samples
        for am2 = 1:number_samples
            sum_joint_probability(am1,am2) = sum(sum(joint_probability{am1,am2}));
        end
    end
%     if round(sum_joint_probability, 4) ~= 1
%         throw('Sum of joint probabilities is not equal to 1.');
%     end
    
    % Calculando a informacao mutua
    MI{1, index} = mutualInformation(size(discretized_noise, 2), number_intervals, marginal_probability, ...
                                     joint_probability, oc);
%     MIaux = zeros(tam,tam);
%     for am1 = 1:tam
%         for am2 = 1:tam
%             for x = 1:number_intervals
%                 for y = 1:number_intervals
%                     fprintf(['# Mutual Information \nOccupancy: ' int2str(oc) ',\t Sample1: ' int2str(am1)  ',\t Sample2: ' int2str(am2) '\n']);
%                     fprintf(['x = ' int2str(x) ',\t y = ' int2str(y) '\n']);
% 
%                     if joint_probability{am1,am2}(x,y) == 0 % necessaria para evitar erros de NaN
%                         continue
%                     end
% 
%                     multprobmarg = marginal_probability(x,am1)*marginal_probability(y,am2);
%                     divprobs = joint_probability{am1,am2}(x,y)/multprobmarg;
%                     MIaux(am1,am2) = MIaux(am1,am2) + joint_probability{am1,am2}(x,y)*log2(divprobs);
%                 end
%             end            
%         end
%     end
%     MI{1, index} = MIaux;
    
    % calculando o crosstalk, que eh a IM contida em todo o processo
    MItotal(index) = mutualInformationCrosstalk(MI{1, index}, number_samples);
%     somamatriztrisup = 0; %soma da matriz triangular superior sem diagonal
%     ini = 2; %variavel da coluna da matriz
%     for i = 1:(tam-1)
%         for j = ini:tam
%             somamatriztrisup = somamatriztrisup + MI{1,index}(i,j);
%         end
%         ini = ini + 1;
%     end
%     somamatriztrisupdiag = 0; %soma da matriz triangular superior com diagonal
%     ini = 1; %variavel da coluna da matriz
%     for i = 1:tam
%         for j = ini:tam
%             somamatriztrisupdiag = somamatriztrisupdiag + MI{1,index}(i,j);
%         end
%         ini = ini + 1;
%     end
%     MItotal(index) = (somamatriztrisup/somamatriztrisupdiag)*100;
end

%Salvando os dados de Informacao Mutua em arquivos .txt
% if semICA
%     dlmwrite('resultadosMI-ICA/MItotal_semICA.txt', MItotal);
%     
%     for oc = 10:10:100
%         ind = (oc/10) + 1;
%         dlmwrite(['resultadosMI-ICA/MIocup' int2str(oc) '_semICA.txt'], MI{1, ind});
%     end
% else
%     dlmwrite('resultadosMI-ICA/MItotal_comICA.txt', MItotal);
%     
%     for oc = 10:10:100
%         ind = (oc/10) + 1;
%         dlmwrite(['resultadosMI-ICA/MIocup' int2str(oc) '_comICA.txt'], MI{1, ind});
%     end
% end

% Plotando os gráficos individualmente
figure
plot(occupancies, MItotal((occupancies/10 + 1), 1), 'Color', [0.6 0.6 0.6], 'Marker', '.', 'MarkerSize', 20, 'MarkerEdgeColor', [0 0.4470 0.7410]);
xlim([0 100]);
title(['Total MI (mPu' int2str(mPu) ')']);
legend('MI before ICA', 'Location', 'northeast');
xlabel('Occupancy (%)');
ylabel('Total Mutual Information (ADC counts)');


% Plotando os graficos antes e depois da ICA
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


