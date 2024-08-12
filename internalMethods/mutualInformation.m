function mutual_information = mutualInformation(number_samples, number_intervals, marginal_probability, ...
                                                joint_probability, occupancy)    
    mutual_information = zeros(number_samples, number_samples);
    for sample1 = 1:number_samples
        for sample2 = 1:number_samples
            for x = 1:number_intervals
                for y = 1:number_intervals
                    fprintf(['# Mutual Information \nOccupancy: ' int2str(occupancy) ...
                             ',\t Sample1: ' int2str(sample1)  ',\t Sample2: ' int2str(sample2) '\n']);

                    fprintf(['x = ' int2str(x) ',\t y = ' int2str(y) '\n']);

                    if joint_probability{sample1, sample2}(x,y) == 0 % necessaria para evitar erros de NaN
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