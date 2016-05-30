%Daniel Balch, Zakari Abner, Tay Johnson
%ECE 449 (Senior Design Project; Adviser: Dr. Dmitrity Garmatyuk)
%2-29-2016
%2-D Positioning Code Deviation Plotter
numSamples = 20e3;
rangesA = [4.731 3.315 3.75 6.12 3.525 2.235 3.12 3.12 1.35 3.06 1.89 1.53 2.415];
rangesB = [2.91 3.315 3.75 3.87 3.525 2.235 3.12 3.12 4.41 3.06 2.25 2.325 2.415];
prePrevChange = 0;
prevChange = 0;
minimaFound = 0;
k = 3;
minimaK = 0;
for i = 1:2
    if (i == 1)
        ranges = rangesA;
    else
        ranges = rangesB;
    end
    while ((k <= (length(ranges))) && (minimaFound == 0))
        prePreviousR = ranges((k-2));
        prevR = ranges((k-1));
        currentR = ranges(k);
        prePrevChange = prevR - prePreviousR;
        prevChange = currentR - prevR;
        cond1 = (prePrevChange <= 0) && (prevChange > 0);
        cond2 = (prePrevChange < 0) && (prevChange >= 0);
        if (cond1 || cond2)
            minimaFound = 1;
            observedY_t = prevR;
            minimaK = (k-1);
        end
        k = k + 1;
    end
    initialR = ranges(1);
    observedX_t = sqrt((initialR^2) - (observedY_t^2));
    observedX_t = abs(observedX_t);
    observedX_p = (sqrt((ranges.^2) - (observedY_t^2)));
    observedX_p = abs(observedX_p);
    for n = 1:(length(observedX_p))
        if (n <= minimaK)
            observedX_p(n) = (-1 * observedX_p(n));
        end
    end
    actualPositions = -36:6:36;
    actualPositions = (actualPositions.*2.54.*0.01);
    deviation =  observedX_p - actualPositions;
    groundDevs = zeros(1,(length(actualPositions)));
    if (i == 1)
        plot(actualPositions,groundDevs);
        hold on;
        plot(actualPositions,deviation,'g','Marker','o');
    else
        plot(actualPositions,deviation,'r','Marker','x');
    end
    if (i == 2)
        legend('Measured Position','Calculated X_{p} (original data)',...
            'Calculated X_{p} (filtered data)','Location','NorthEast');
        grid on;
        xlabel('Actual Horizontal Position (meters; negative means left of center)', 'FontSize', 14);
        ylabel('Deviation (meters)', 'FontSize', 14);
        title('Deviations from Actual Position using Original and Filtered Data', 'FontSize', 18);
    end
end