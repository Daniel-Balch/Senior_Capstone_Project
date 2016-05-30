%Zak Abner, Tay Johnson, Daniel Balch
%ECE 449 A (Professor Dmitriy Garmatyuk)
%Best-fit Navigation Algorithm for 2-29 Experiment
%This algorithm assumes that the length of rangesA is the same as the
%length of ranges B, and that there are at least 3 data points in each.
%By convention, the first data point will have an x-coordinate of 0.
%(X_t is no longer preset to 0.)
rangesA = [4.731 3.315 3.75 6.12 3.525 2.235 3.12 3.12 1.35 3.06 1.89 1.53 2.415];
rangesB = [6.351 7.68 3.99 6.21 3.57 7.38 4.155 3.12 4.65 3.06 5.115 3.915 4.74];
posInterval = 6; % acutal position interval in inches
centerIndex = 7; % index of actual center position (closest to radar)
actualPositions = 0:posInterval:(((length(rangesA))-1)*posInterval);
actualPositions = (actualPositions.*2.54.*0.01); % actual positions in meters
actualY_t = 118*2.54*0.01;
actualX_t = actualPositions(centerIndex);
actualRanges = sqrt(((actualPositions-actualX_t).^2)+(actualY_t^2));
rngDevsA = (rangesA - actualRanges);
rngDevsB = (rangesB - actualRanges);
rLength = length(rangesA);
observedX_pA = zeros((rLength - 2),rLength);
observedX_pB = zeros((rLength - 2),rLength);
posDevsA = zeros((rLength - 2),rLength);
posDevsB = zeros((rLength - 2),rLength);
observedX_tA = zeros(1,(rLength - 2));
observedX_tB = zeros(1,(rLength - 2));
x_tDevsA = zeros(1,(rLength - 2));
x_tDevsB = zeros(1,(rLength - 2));
observedY_tA = zeros(1,(rLength - 2));
observedY_tB = zeros(1,(rLength - 2));
y_tDevsA = zeros(1,(rLength - 2));
y_tDevsB = zeros(1,(rLength - 2));
parabCoeffsA = zeros((rLength-2),3);
parabCoeffsB = zeros((rLength-2),3);
invalPreMarkersA = zeros((rLength - 2),rLength);
invalPreMarkersB = zeros((rLength - 2),rLength);
parabola = @(coefs,x)((coefs(1).*(x.^2))+(coefs(2).*x)+coefs(3));
initCoefs = zeros(1,3);
normResidualsA = zeros(1,(rLength - 2));
normResidualsB = zeros(1,(rLength - 2));
for k = 3:(rLength)
    k1 = (k-2);
    for u = 1:2
        if (u == 1)
            currentRngs = rangesA;
        else
            currentRngs = rangesB;
        end
        currentIndices = (1:k);
        currentY = currentRngs(1:k);
        [coefs,currentNormResidual] = lsqcurvefit(parabola,initCoefs,...
            currentIndices,currentY);
        currentMinimaIndex = (-1*(coefs(2)/(2*coefs(1))));
        currentY_t = parabola(coefs,currentMinimaIndex);
        currentX_t = abs(sqrt((((currentRngs(1))^2)-(currentY_t^2))));
        if (currentMinimaIndex < 1)
            currentX_t = (-1*currentX_t);
        end
        currentObsX_p = zeros(1,length(currentRngs));
        currentInvalPreMarkers = zeros(1,length(currentRngs));
        for m = 1:(length(currentRngs))
            preX_p = (((currentRngs(m))^2)-(currentY_t^2));
            if (m >= currentMinimaIndex)
                currentObsX_p(m) = (currentX_t+(abs(sqrt(preX_p))));
            else
                currentObsX_p(m) = (currentX_t-(abs(sqrt(preX_p))));
            end
            if (preX_p < 0)
                currentInvalPreMarkers(m) = 1;
            end
        end
        if (u == 1)
            observedX_pA(k1,:) = currentObsX_p;
            posDevsA(k1,:) = (currentObsX_p - actualPositions);
            observedX_tA(k1) = currentX_t;
            x_tDevsA(k1) = (currentX_t - actualX_t);
            observedY_tA(k1) = currentY_t;
            y_tDevsA(k1) = (currentY_t - actualY_t);
            parabCoeffsA(k1,:) = coefs;
            invalPreMarkersA(k1,:) = currentInvalPreMarkers;
            normResidualsA(k1) = currentNormResidual;
        else
            observedX_pB(k1,:) = currentObsX_p;
            posDevsB(k1,:) = (currentObsX_p - actualPositions);
            observedX_tB(k1) = currentX_t;
            x_tDevsB(k1) = (currentX_t - actualX_t);
            observedY_tB(k1) = currentY_t;
            y_tDevsB(k1) = (currentY_t - actualY_t);
            parabCoeffsB(k1,:) = coefs;
            invalPreMarkersB(k1,:) = currentInvalPreMarkers;
            normResidualsB(k1) = currentNormResidual;
        end
    end
end
clear currentRngs currentX currentY coefs currentNormResidual currentMinimaIndex;
clear currentX_t currentY_t currentObsX_p currentInvalPreMarkers preX_p;
fgCount = 1;
zeroDev = zeros(1,100);
actPosMargin = (((max(actualPositions))-(min(actualPositions)))*0.05);
minActPos = ((min(actualPositions))-actPosMargin);
minParabActPos = ((actualPositions(3))-actPosMargin);
maxActPos = ((max(actualPositions))+actPosMargin);
topMarginFactor = 6;
parabolaX = linspace(minActPos,maxActPos,500);
parabolaIndices = (((1/(posInterval*2.54*0.01)).*(parabolaX))+1);
for v = 1:2
    if (v == 1)
        currentPosDevs = posDevsA;
        currentX_tDevs = x_tDevsA;
        currentY_tDevs = y_tDevsA;
        currentPreMarkers = invalPreMarkersA;
        currentRngDevs = rngDevsA;
        currentRngs = rangesA;
        currentParabCoeffs = parabCoeffsA;
        typeText = 'original data';
    else
        currentPosDevs = posDevsB;
        currentX_tDevs = x_tDevsB;
        currentY_tDevs = y_tDevsB;
        currentPreMarkers = invalPreMarkersB;
        currentRngDevs = rngDevsB;
        currentRngs = rangesB;
        currentParabCoeffs = parabCoeffsB;
        typeText = 'filtered data';
    end
    figure(fgCount);
    fgCount = (fgCount+1);
    scatter(actualPositions(3:(rLength)),currentX_tDevs,'g');
    hold on;
    plot(actualPositions(3:(rLength)),zeroDev(1:(rLength-2)),'b');
    xlabel('Horizontal Position (meters)','FontSize',14);
    ylabel('X_{T} Deviation (meters)','FontSize',14);
    titleText = sprintf('Best-Fit Curve X_{T} Deviation vs Position (%s)',typeText);
    title(titleText,'FontSize',18);
    %fprintf('\nX_t Legend, Figure %i, v = %i',(fgCount-1),v);
    legendText = sprintf('X_{T} Deviation (%s)',typeText);
    legend(legendText,'Actual X_{T} (no deviation)','Location','NorthEast');
    maxYCoord = max([currentX_tDevs 0]);
    minYCoord = min([currentX_tDevs 0]);
    yMargin = ((maxYCoord-minYCoord)*0.1);
    minYCoord = (minYCoord-yMargin);
    maxYCoord = (maxYCoord+(yMargin*topMarginFactor));
    axis([minParabActPos,maxActPos,minYCoord,maxYCoord]);
    grid on;
    figure(fgCount);
    fgCount = (fgCount+1);
    scatter(actualPositions(3:(rLength)),currentY_tDevs,'g');
    hold on;
    plot(actualPositions(3:(rLength)),zeroDev(1:(rLength-2)),'b');
    xlabel('Horizontal Position (meters)','FontSize',14);
    ylabel('Y_{T} Deviation (meters)','FontSize',14);
    titleText = sprintf('Best-Fit Curve Y_{T} Deviation vs Position (%s)',typeText);
    title(titleText,'FontSize',18);
    %fprintf('\nY_t Legend, Figure %i, v = %i',(fgCount-1),v);
    legendText = sprintf('Y_{T} Deviation (%s)',typeText);
    legend(legendText,'Actual Y_{T} (no deviation)','Location','NorthEast');
    maxYCoord = max([currentY_tDevs 0]);
    minYCoord = min([currentY_tDevs 0]);
    yMargin = ((maxYCoord-minYCoord)*0.1);
    minYCoord = (minYCoord-yMargin);
    maxYCoord = (maxYCoord+(yMargin*topMarginFactor));
    axis([minParabActPos,maxActPos,minYCoord,maxYCoord]);
    grid on;
    figure(fgCount);
    fgCount = (fgCount+1);
    scatter(actualPositions,currentRngDevs,'g');
    hold on;
    plot(actualPositions,zeroDev(1:(length(actualPositions))),'b');
    xlabel('Horizontal Position (meters)','FontSize',14);
    ylabel('Range Deviation (meters)','FontSize',14);
    titleText = sprintf('Range Deviations by Position (%s)',typeText);
    title(titleText,'FontSize',18);
    %fprintf('\nRange Dev. Legend, Figure %i, v = %i',(fgCount-1),v);
    legendText = sprintf('Range Deviations (%s)',typeText);
    legend(legendText,'Actual Ranges (no deviation)','Location','NorthEast');
    maxYCoord = max([currentRngDevs 0]);
    minYCoord = min([currentRngDevs 0]);
    yMargin = ((maxYCoord-minYCoord)*0.1);
    minYCoord = (minYCoord-yMargin);
    maxYCoord = (maxYCoord+(yMargin*topMarginFactor));
    axis([minActPos,maxActPos,minYCoord,maxYCoord]);
    grid on;
    for n = 1:(rLength-2)
        curveNum = (n+2);
        preMarkerSeg = currentPreMarkers(n,:);
        posDevSeg = currentPosDevs(n,:);
        numInvMarkers = sum(preMarkerSeg);
        parabCoeffsSeg = currentParabCoeffs(n,:);
        parabolaY = parabola(parabCoeffsSeg,parabolaIndices);
        if (numInvMarkers > 0)
            invMarkerIndices = ones(1,numInvMarkers);
            invMarkerCount = 1;
            for z = 1:(length(preMarkerSeg))
                if (preMarkerSeg(z) == 1)
                    invMarkerIndices(invMarkerCount) = z;
                    invMarkerCount = (invMarkerCount+1);
                end
            end
            invMarkerX = actualPositions(invMarkerIndices);
            invMarkerY = posDevSeg(invMarkerIndices);
            invMarkerR = currentRngs(invMarkerIndices);
        end
        figure(fgCount);
        fgCount = (fgCount+1);
        scatter(actualPositions,posDevSeg,'g');
        hold on;
        plot(actualPositions,zeroDev(1:(length(actualPositions))),'b');
        if (numInvMarkers > 0)
            scatter(invMarkerX,invMarkerY,'r','filled');
        end
        xlabel('Actual Horizontal Position (meters)','FontSize',12);
        ylabel('X_{P} Deviation (meters)','FontSize',12);
        if (curveNum < (length(actualPositions)))
            titleText = sprintf('X_{P} Deviations for Best-Fit Curve Generated at x = %.2f meters (%s)',...
                actualPositions(curveNum),typeText);
        else
            titleText = sprintf('X_{P} Deviations for Final Best-Fit Curve (%s)',...
                typeText);
        end
        title(titleText,'FontSize',16);
        %fprintf('\n\nX_p Legend, Parabola %i, Figure %i, v = %i',...
            %n,(fgCount-1),v);
        legendText = sprintf('X_{P} Deviation (%s)',typeText);
        if (numInvMarkers > 0)
            legend(legendText,'Actual X_{P} (no deviation)',...
                'Invalid/Impossible Data Points','Location','NorthEast');
        else
            legend(legendText,'Actual X_{P} (no deviation)','Location',...
                'NorthEast');
        end
        maxYCoord = max([posDevSeg 0]);
        minYCoord = min([posDevSeg 0]);
        yMargin = ((maxYCoord-minYCoord)*0.1);
        minYCoord = (minYCoord-yMargin);
        maxYCoord = (maxYCoord+(yMargin*topMarginFactor));
        axis([minActPos,maxActPos,minYCoord,maxYCoord]);
        grid on;
        figure(fgCount);
        fgCount = (fgCount+1);
        scatter(actualPositions,currentRngs,'b');
        hold on;
        plot(parabolaX,parabolaY,'g');
        if (numInvMarkers > 0)
            scatter(invMarkerX,invMarkerR,'r','filled');
        end
        xlabel('Actual Horizontal Position (meters)','FontSize',12);
        ylabel('Measured Range (meters)','FontSize',12);
        if (curveNum < (length(actualPositions)))
            titleText = sprintf('Ranges with Best-Fit Curve Generated at x = %.2f meters (%s)',...
                actualPositions(curveNum),typeText);
        else
            titleText = sprintf('Ranges with Final Best-Fit Curve (%s)',...
                typeText);
        end
        title(titleText,'FontSize',16);
        %fprintf('\n\nX_p Legend, Parabola %i, Figure %i, v = %i',...
            %n,(fgCount-1),v);
        legendText = sprintf('Radar Range in meters (%s)',typeText);
        if (numInvMarkers > 0)
            legend(legendText,'Best-Fit Curve of Ranges',...
                'Invalid/Impossible Data Points','Location','NorthEast');
        else
            legend(legendText,'Best-Fit Curve of Ranges','Location',...
                'NorthEast');
        end
        maxYCoord = max(currentRngs);
        minYCoord = min(currentRngs);
        yMargin = ((maxYCoord-minYCoord)*0.1);
        minYCoord = (minYCoord-yMargin);
        maxYCoord = (maxYCoord+(yMargin*topMarginFactor));
        axis([minActPos,maxActPos,minYCoord,maxYCoord]);
        grid on;
    end
end
% post matrix creation:
%   -for each parabola, make array of invalid-marker x-coordinates
%   -use loop to separate ranges from positions;
%   -make array of invalid y-coordinates for each of above points
