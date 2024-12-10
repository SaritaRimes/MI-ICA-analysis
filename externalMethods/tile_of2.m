% implementa OF2 com 3 restrições
% bernardo smp + sarita mr

function [weights] = tile_of2(noise, optimization)
    % Defining the pulse and the pulse derivative
    pulse = [0.0000 0.0172 0.4524 1.0000 0.5633 0.1493 0.0424];
    d_pulse = [0.00004019 0.00333578 0.03108120 0.00000000 -0.02434490 -0.00800683 -0.00243344];
    
    % Checking if the choice is optimize or not
    if optimization == 0 || ~optimization
	    covariance = eye(7,7);
    elseif optimization == 1 || optimization
	    covariance = cov(noise);
    else
        throw("Options for ""optimize"" are 0 (False) or 1 (True).");
    end

    % Making the B matrix
    b = [covariance; pulse; d_pulse; ones(1,7)];
    b = [b [pulse'; zeros(3,1)] [d_pulse'; zeros(3,1)] [ones(7,1); zeros(3,1)]];
    
    % Making the C matrix
    c = [zeros(7,1);1;zeros(2,1)];
    
    % Calculate the weights
    aux = b\c;
    weights = aux(1:7)';
end
