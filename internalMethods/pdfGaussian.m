function pdf_gaussian = pdfGaussian(mean, covariance, r, s, A)
    if all(size(A) == 1) % A is a single number
        interestSignal = A*s;
    elseif (isrow(A) || iscolumn(A)) && ismember(7, size(A))
        interestSignal = A;
    else
        throw(['The array of amplitudes does not match the size required. ' ...
               'It must be a row or columns array, or a matrix with ' ...
               'one of the two dimensions equal to 7.']);
    end

    fprintf('Estimating Gaussian PDF...');

    noise = r - interestSignal;

    pdf_gaussian = (1/(((2*pi)^3.5)*sqrt(det(covariance)))) ...
                  *exp((-(noise - mean)*inv(covariance) ...
                       *(noise - mean)')/2);
end