function chi2 = chi2EfficiencyEstimation(amplitude, r, s, pedestal)
    noise_temporary = r - (amplitude*s + pedestal);

    chi2 = sqrt(sum((noise_temporary.^2)./7));
end