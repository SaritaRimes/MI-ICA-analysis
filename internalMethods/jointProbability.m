function joint_prob = jointProbability(noise, number_intervals, occupancy)
    fprintf(['# Joint Probability \nOccupancy: ' int2str(occupancy) '\n']);

    discretized_noise = noise; % serao as amostras discretizadas

    number_events = size(discretized_noise, 1); % number of events in data
    number_samples = size(discretized_noise, 2); % number of samples in data

    joint_prob = cell(number_samples, number_samples);
    for sample1 = 1:number_samples
        for sample2 = 1:number_samples
            joint_prob_aux = zeros(number_intervals, number_intervals);

            for x = 1:number_intervals
                for y = 1:number_intervals
                    for i = 1:number_events
                        if (discretized_noise(i, sample1) == x && discretized_noise(i, sample2) == y)
                            joint_prob_aux(x,y) = joint_prob_aux(x,y) + 1;
                        end
                    end
                end
            end

            joint_prob{sample1, sample2} = joint_prob_aux/number_events;
        end
    end
end