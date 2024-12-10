function [probabilities, chi2, stdErrorCutsProb, stdErrorCutsChi2] ...
            = analysisDataQuality(amplitudes, receivedSignal, ...
                                  meanLognormal, meanGaussian, ...
                                  covarianceLognormal, covarianceGaussian, ...
                                  pedestal, occupancy)

    addpath('externalMethods', 'estimationMethods');

    % Defining the amplitudes
    amplitudeOF2 = amplitudes{1};
    amplitudeGauss = amplitudes{2};
    amplitudeCOFAll = amplitudes{3};
    amplitudeCOF = amplitudeCOFAll(:, 4);
    amplitudeLogn = amplitudes{4};
    amplitudeTrue = amplitudes{5};

    eventsNumber = length(amplitudeTrue);

    % Defining some constants
    cutsNumber = 100000;
    s = [0 .0172 .4524 1 .5633 .1493 .0424]; % calorimeter reference pulse

    % Initializing some structures
    pdfOF2 = zeros(eventsNumber, 1);
    pdfGauss = zeros(eventsNumber,1);
    pdfCOF = zeros(eventsNumber,1);
    pdfLogn = zeros(eventsNumber,1);
    pdfLognGauss = zeros(eventsNumber,1);
    chi2OF2 = zeros(eventsNumber, 1);
    chi2Gauss = zeros(eventsNumber, 1);
    chi2COF = zeros(eventsNumber, 1);
    chi2Logn = zeros(eventsNumber, 1);

    % Estimating the values of probabilities and chi2
    H = load('matrizCOF.txt');
    for i = 1:eventsNumber

        fprintf(['Probabilities and chi2, Iteration ' int2str(i) '/' ...
                 int2str(eventsNumber) '\n']);
      
        % PDF OF2
        pdfOF2(i) = pdfgaussian(meanGaussian, covarianceGaussian, ...
                                receivedSignal(i,:), s, amplitudeOF2(i));

        % PDF Gaussian
        pdfGauss(i) = pdfgaussian(meanGaussian, covarianceGaussian, ...
                                  receivedSignal(i,:), s, ...
                                  amplitudeGauss(i));

        % PDF COF
        pdfCOF(i) = pdfgaussian(meanGaussian, covarianceGaussian, ...
                                receivedSignal(i,:), s, ...
                                amplitudeCOFAll(i,:)*H');

        % PDF Lognormal
        pdfLognGauss(i) = pdfgaussian(meanGaussian, covarianceGaussian, ...
                                      receivedSignal(i,:), s, ...
                                      amplitudeLogn(i));
        pdfLogn(i) = pdflognormal(meanLognormal, covarianceLognormal, ...
                                  receivedSignal(i,:), s, amplitudeLogn(i));

        % Chi2 OF2
        amp = amplitudeOF2(i)*s + pedestal;
        chi2OF2(i) = sqrt(sum(((receivedSignal(i,:) - amp).^2)./7));

        % Chi2 Gaussian
        amp = amplitudeGauss(i)*s + pedestal;
        chi2Gauss(i) = sqrt(sum(((receivedSignal(i,:) - amp).^2)./7));

        % Chi2 COF
        amp = amplitudeCOF(i)*s + pedestal;
        chi2COF(i) = sqrt(sum(((receivedSignal(i,:) - amp).^2)./7));
    
        % Chi2 Lognormal
        amp = amplitudeLogn(i)*s + pedestal;
        chi2Logn(i) = sqrt(sum(((receivedSignal(i,:) - amp).^2)./7));

    end

    % Calculating the step of each method
    % OF2
    stepOF2Prob = (max(pdfOF2) - min(pdfOF2))/cutsNumber;
    stepOF2Chi2 = (max(chi2OF2) - min(chi2OF2))/cutsNumber;
    % MLE Gaussiano
    stepGaussProb = (max(pdfGauss) - min(pdfGauss))/cutsNumber;
    stepGaussChi2 = (max(chi2Gauss) - min(chi2Gauss))/cutsNumber;
    % COF
    stepCOFProb = (max(pdfCOF) - min(pdfCOF))/cutsNumber;
    stepCOFChi2 = (max(chi2COF) - min(chi2COF))/cutsNumber;
    % MLE Lognormal
    switch occupancy
        case 10
            stepLognProb = (max(pdfLogn) - min(pdfLogn))/(cutsNumber); % 10%
        case 30
            stepLognProb = (max(pdfLogn) - min(pdfLogn))/(cutsNumber); % 30%
        case 50
            stepLognProb = (max(pdfLogn) - min(pdfLogn))/(cutsNumber); % 50%
        case 80
            stepLognProb = (max(pdfLogn) - min(pdfLogn))/(cutsNumber); % 80%
        otherwise
            throw(MException('Enter a valid occupancy. The options are: 10, 30, 50 or 80.'));
    end
    stepLognChi2 = (max(chi2Logn) - min(chi2Logn))/cutsNumber;

    % Estimating the error
    erroGauss = amplitudeGauss - amplitudeTrue;
    erroOF2 = amplitudeOF2 - amplitudeTrue;
    erroCOF = amplitudeCOF - amplitudeTrue;
    erroLogn = amplitudeLogn - amplitudeTrue;

    % Defining some structures
    stdErrorCutsOF2Prob = zeros(cutsNumber, 1);
    stdErrorCutsGaussProb = zeros(cutsNumber, 1);
    stdErrorCutsCOFProb = zeros(cutsNumber, 1);
    stdErrorCutsLognProb = zeros(cutsNumber, 1);
    stdErrorCutsOF2Chi2 = zeros(cutsNumber, 1);
    stdErrorCutsGaussChi2 = zeros(cutsNumber, 1);
    stdErrorCutsCOFChi2 = zeros(cutsNumber, 1);
    stdErrorCutsLognChi2 = zeros(cutsNumber, 1);

    % Applying the cuts and recalculating the std of chi2 and probability
    j = cutsNumber;
    for i = 0:1:cutsNumber

        fprintf(['Cuts, Iteration: ' int2str(i) '/' int2str(cutsNumber) '\n']);
    
        % Probability OF2
        indOF2 = pdfOF2 > stepOF2Prob*(i);
        stdErrorCutsOF2Prob(i + 1) = std(erroOF2(indOF2));

        % Probability Gaussiano
        indGaus = pdfGauss > stepGaussProb*(i);
        stdErrorCutsGaussProb(i + 1) = std(erroGauss(indGaus));
    
        % Probability COF
        indCOF = pdfCOF > stepCOFProb*(i);
        stdErrorCutsCOFProb(i + 1) = std(erroCOF(indCOF));
    
        % Probability Lognormal
        indLogn = pdfLogn > stepLognProb*(i);
        stdErrorCutsLognProb(i + 1) = std(erroLogn(indLogn));
    
        % Chi2 OF2
        indOF2 = chi2OF2 < stepOF2Chi2*(j);
        stdErrorCutsOF2Chi2(i + 1) = std(erroOF2(indOF2));
    
        % Chi2 Gaussiano
        indGaus = chi2Gauss < stepGaussChi2*(j);
        stdErrorCutsGaussChi2(i + 1) = std(erroGauss(indGaus));
    
        % Chi2 COF
        indCOF = chi2COF < stepCOFChi2*(j);
        stdErrorCutsCOFChi2(i + 1) = std(erroCOF(indCOF));
    
        % Chi2 Lognormal
        indLogn = chi2Logn < stepLognChi2*(j);
        stdErrorCutsLognChi2(i + 1) = std(erroLogn(indLogn));

        j = j - 1;
    
    end

    probabilities = {pdfGauss, pdfOF2, pdfCOF, pdfLogn, pdfLognGauss};
    chi2 = {chi2Gauss, chi2OF2, chi2COF, chi2Logn};
    stdErrorCutsProb = {stdErrorCutsGaussProb, stdErrorCutsOF2Prob, ...
                        stdErrorCutsCOFProb, stdErrorCutsLognProb};
    stdErrorCutsChi2 = {stdErrorCutsGaussChi2, stdErrorCutsOF2Chi2, ...
                        stdErrorCutsCOFChi2, stdErrorCutsLognChi2};

end