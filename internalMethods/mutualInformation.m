function mutual_information = mutualInformation(number_samples, number_intervals, marginal_probability, ...
                                                joint_probability, occupancy, fold)
    fprintf(['# Mutual Information \nOccupancy: ' int2str(occupancy) ', \tFold: ' int2str(fold) '\n']);
    
    mutual_information = zeros(number_samples, number_samples);
    for sample1 = 1:number_samples
        for sample2 = 1:number_samples
            for x = 1:number_intervals
                for y = 1:number_intervals
                    if joint_probability{sample1, sample2}(x,y) == 0 % necessary to avoid NaN errors
                        continue
                    end

                    multprobmarg = marginal_probability(x, sample1)*marginal_probability(y, sample2);
                    divprobs = joint_probability{sample1, sample2}(x,y)/multprobmarg;
                    mutual_information(sample1, sample2) = mutual_information(sample1, sample2) + joint_probability{sample1, sample2}(x,y)*log2(divprobs);
                end
            end            
        end
    end
end