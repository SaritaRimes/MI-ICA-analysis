function checkProbabilities(probability_type, probability)
    if probability_type == "marginal"
        sum_probability = sum(probability);

        if round(sum_probability, 4) ~= 1
            throw('Sum of marginal probabilities is not equal to 1.');
        end
    elseif probability_type == "joint"
        number_samples = size(probability, 1);
        sum_probability = zeros(number_samples, number_samples);

        for sample1 = 1:number_samples
            for sample2 = 1:number_samples
                sum_probability(sample1, sample2) = sum(sum(probability{sample1, sample2}));
            end
        end
        if round(sum_probability, 4) ~= 1
            throw('Sum of joint probabilities is not equal to 1.');
        end
    else
        throw('Probability type is not valid. The options are ''marginal'' or ''joint''.');
    end
end