function amplitude_ica = amplitudeIca(r, s, mean_noise_training, ica_matrix, amplitude_gauss, ...
                                      number_dimensions, marginal_probability, spline_hist)
    
    maximum_probability = -1;
    for amplitude_auxiliary = (amplitude_gauss - 100):1:(amplitude_gauss + 100)
        
        noise_auxiliary_temporary = r - amplitude_auxiliary*s;
        noise_auxiliary = (noise_auxiliary_temporary - mean_noise_training)*ica_matrix';
        
        for j = 1:number_dimensions
            marginal_probability(j) = ppval(spline_hist(j), noise_auxiliary(j));
        end
        
        probability_ica = prod(marginal_probability);
        
        if (probability_ica > maximum_probability)
            amplitude = amplitude_auxiliary;
            maximum_probability = probability_ica;
        end
    end
    
    amplitude_ica = amplitude;
end