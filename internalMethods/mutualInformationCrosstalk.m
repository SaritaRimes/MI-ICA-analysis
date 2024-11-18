function mutual_information_crosstalk = mutualInformationCrosstalk(mutual_information, number_samples)
    upper_triangular_matrix_sum = 0; % sum of superior triangular matrix without its diagonal
    ini = 2; % variable of column of matrix
    for i = 1:(number_samples-1)
        for j = ini:number_samples
            upper_triangular_matrix_sum = upper_triangular_matrix_sum + mutual_information(i,j);
        end

        ini = ini + 1;
    end

    lower_triangular_matrix_sum = 0; % sum of lower triangular matrix with its diagonal
    ini = 1; % variable of column of matrix
    for i = 1:number_samples
        for j = ini:number_samples
            lower_triangular_matrix_sum = lower_triangular_matrix_sum + mutual_information(i,j);
        end

        ini = ini + 1;
    end

    mutual_information_crosstalk = (upper_triangular_matrix_sum/lower_triangular_matrix_sum)*100;
end