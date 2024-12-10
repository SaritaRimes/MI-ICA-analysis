function pdfLognormal = pdflognormal(mean, covariance, r, s, A)
    noise = r - A*s;
    if any(noise <= 0)
        pdfLognormal = 0;
    else
        numerator = exp(-.5*((log(noise) - mean)*inv(covariance) ...
                    *(log(noise) - mean)'));
        denominator = prod(noise)*sqrt(det(covariance))*(2*pi)^(3.5);
    
        pdfLognormal = numerator / denominator;
    end
end