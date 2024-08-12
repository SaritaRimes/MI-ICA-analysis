function [discretized_samples, number_samples_interval] = discretizeSamples(samples, number_intervals, occupancy)
    discretized_samples = samples; % discretized samples

    number_samples = size(discretized_samples, 2); % number of samples in data
    step = (1-(-1))/number_intervals;

    number_samples_interval = -1*ones(number_intervals, number_samples);
    for samp = 1:number_samples % iterates in samples
        inte = 1; % to count the intervals

        for i = -1:step:(1 - step)
            fprintf(['# Discretization \nOccupancy: ' int2str(occupancy) ...
                     ',\t Sample: ' int2str(samp)  ',\t Interval: ' num2str(i) '\n']);

            if i == 1 - step
                k = find((samples(:, samp) >= i) & (samples(:, samp) <= (i + step)));
                discretized_samples(k, samp) = inte;
            else
                k = find((samples(:, samp) >= i) & (samples(:, samp) < (i + step)));
                discretized_samples(k, samp) = inte;
            end

            number_samples_interval(inte, samp) = size(k, 1);
            inte = inte + 1;
            clear k
        end
    end
end