function [amplitude_ica, probability_ica] = IcaEstimator(r, s, noise_ica, mean_noise_training, ica_matrix, ...
                                                         number_dimensions, number_bins, amplitude_gauss, ...
                                                         occupancy, iteration)

    % Normalizing the histograms of noiseand finding their x coordinates
    hist_probabilities = -1*ones(number_bins, size(noise_ica, 2));
    hist_bins = zeros(number_bins + 1, size(noise_ica,2));
    hist_coordinate_x = zeros(size(hist_bins, 1) - 1, number_dimensions);
    ica_probabilities_hist = gobjects(1, number_dimensions);
    for i = 1:number_dimensions
        % Plotting the normalized histograms
        ica_probabilities_hist(i) = figure;
        h = histogram(noise_ica(:, i), number_bins, 'Normalization', 'probability', ...
                      'FaceColor', '#DF65F8', 'EdgeColor', '#DF65F8');

        % Finding the bin edges and the values
        hist_probabilities(:, i) = h.Values;
        hist_bins(:, i) = h.BinEdges;

        % Finding the x coordinates
        for j = 1:size(hist_bins, 1) - 1
            hist_coordinate_x(j, i) = (hist_bins(j, i) + hist_bins(j + 1, i))/2;
        end

        % Using splines to interpolate
        spline_hist(i) = spline(hist_coordinate_x(:, i), hist_probabilities(:, i));

        % Plotting the interpolations
        hold on
        plot(hist_coordinate_x(:, i), ppval(hist_coordinate_x(:, i), spline_hist(i)), ...
             'Color', 'k', 'LineWidth', 1);
        legend({['variable ' int2str(i)], ['interpolation ' int2str(i)]}, 'Location', 'best');
        hold off
    end

    % Estimating the amplitude using MLE + ICA method
    amplitude_ica = amplitude_gauss;
    number_events = size(amplitude_ica, 1);
    marginal_probability = zeros(1, number_dimensions);
    probability_ica = ones(number_events, 1);
    for i = 1:number_events
        fprintf(['MLE + ICA: ' ...
                 'Processing event ' int2str(i) '/' int2str(number_events) '...' ...
                 '\nOccupancy ' int2str(occupancy) ', iteration ' int2str(iteration) '\n']);

        amplitude_ica(i) = amplitudeIca(r(i,:), s, mean_noise_training, ica_matrix, amplitude_gauss(i), ...
                                        number_dimensions, marginal_probability, spline_hist);
        probability_ica(i) = pdfIca(r(i,:), s, mean_noise_training, ica_matrix, amplitude_gauss(i), ...
                                    number_dimensions, marginal_probability, spline_hist);
    end
end