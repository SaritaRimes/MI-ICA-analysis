function [discretized_samples, number_samples_interval] = discretizeSamples(noise, number_intervals, occupancy)
    discretized_samples = noise; % discretized samples

    number_samples = size(discretized_samples, 2); % number of samples in data
    step = (1-(-1))/number_intervals;

    number_samples_interval = -1*ones(number_intervals, number_samples);
    for sample = 1:number_samples
        inte = 1; % to count the intervals

        for i = -1:step:(1 - step)
            fprintf(['# Discretization \nOccupancy: ' int2str(occupancy) ...
                     ',\t Sample: ' int2str(sample)  ',\t Interval: ' num2str(i) '\n']);

            if i == 1 - step
                k = find((noise(:, sample) >= i) & (noise(:, sample) <= (i + step)));
                discretized_samples(k, sample) = inte;
            else
                k = find((noise(:, sample) >= i) & (noise(:, sample) < (i + step)));
                discretized_samples(k, sample) = inte;
            end

            number_samples_interval(inte, sample) = size(k, 1);
            inte = inte + 1;
            clear k
        end
    end
end