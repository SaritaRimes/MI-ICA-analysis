function x = golden_section(objectiveFunction, xlower, xupper, tolerance)
    % Calculating the golden ratio
    proportion = (sqrt(5) - 1) / 2;

    % Checking value of the function in the bounders
    flower = objectiveFunction(xlower);
    fupper = objectiveFunction(xupper);
    while flower == 0 || fupper == 0
        if flower == 0 && fupper == 0
            xlower = xlower + 1;
            xupper = xupper - 1;
            flower = objectiveFunction(xlower);
            fupper = objectiveFunction(xupper);
        elseif flower == 0
            xlower = xlower + 1;
            flower = objectiveFunction(xlower);
        elseif fupper == 0
            xupper = xupper - 1;
            fupper = objectiveFunction(xupper);
        end

        if xupper <= xlower
            x = 0;
            return
        end
    end

    % Calculating the internal points
    x1 = xlower + proportion * (xupper - xlower);
    x2 = xupper - proportion * (xupper - xlower);

    % Estimating value of the function at the internal points
    fx1 = objectiveFunction(x1);
    fx2 = objectiveFunction(x2);

    % Extimating the parameter of interest
    while xupper - xlower >= tolerance
        if fx1 > fx2
            xlower = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = xlower + proportion * (xupper - xlower);
            fx1 = objectiveFunction(x1);
        else
            xupper = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = xupper - proportion * (xupper - xlower);
            fx2 = objectiveFunction(x2);
        end
    end

    % Estimating the final x
%     flower = objectiveFunction(xlower);
%     fupper = objectiveFunction(xupper);
%     fmid = objectiveFunction((xupper + xlower) / 2);
%     if flower > fupper && flower > fmid
%         x = xlower;
%     elseif fupper > flower && fupper > fmid
%         x = xupper;
%     else
%         x = (xupper + xlower) / 2;
    x = (xupper + xlower) / 2;
end