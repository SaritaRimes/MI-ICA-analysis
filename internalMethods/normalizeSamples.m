function normalized_samples = normalizeSamples(samples)
    number_events = size(samples, 1);
    number_samples = size(samples, 2);

    normalized_samples = zeros(number_events, number_samples);
    for i = 1:number_samples
        nume = 2*(samples(:,i) - min(samples(:,i)));
        deno = max(samples(:,i)) - min(samples(:,i));
        normalized_samples(:,i) = (nume/deno) - 1;
    end
end