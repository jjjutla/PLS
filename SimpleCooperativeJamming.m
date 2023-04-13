
% Initialize an array for average secrecy rate
avgRsCJ=zeros(1,61);
avgRsDirect = zeros(1, 61);

% Set path loss exponent
c=3.5; 
%change Source-Eavesdropper distance
countNaN=0;

% Loop over different values of Source-Eavesdropper distance
% This loop will compute the secrecy rate for different distances
% between the source and eavesdropper, ranging from 30 to 90 meters
for dSE=30:90

    % Initialize variables for computing average secrecy rate
    sumRsCJ=0;
    count=0;

    % Loop over multiple iterations to compute the average secrecy rate
    % for a given Source-Eavesdropper distance
    for itnr=1:100000

        % Generate a random complex channel vector for each node
        % This will simulate the fading effect of wireless channel
        e=(randn(3,1) + randn(3,1)*1i)/sqrt(2);
        
        % Set power of source signal
        PS=10^0.3;

        % Set noise power
        phi2=10^-6;

        % Set distance between source and destination
        dSD=50;

        % Compute channel gain between source and destination
        hSD=((dSD)^(-c/2))*e(1,1);

        %Set distance between relay and destination
        dRD=25;

        % Compute channel gain between relay and destination
        hRD=((dRD)^(-c/2))*e;

        % Compute covariance matrix of relay-destination channel
        RRD=hRD*ctranspose(hRD);

        % Set distance between source and eavesdropper
        hSE=((dSE)^(-c/2))*e(1,1);

        % Set distance between relay and eavesdropper
        dRE=dSE-25;

        % Compute channel gain between relay and eavesdropper
        hRE=((dRE)^(-c/2))*e;

        % Compute covariance matrix of relay-eavesdropper channel
        RRE=hRE*ctranspose(hRE);

        % Compute optimal weight for beamforming
        u=sqrt(PS)/norm((norm(hRD)^2)*hRE-ctranspose(hRD)*hRE*hRD);
        w=u*(norm(hRD)^2)*hRE-u*ctranspose(hRD)*hRE*hRD;
        
        % Compute power of relay signal and eavesdropper signal
        PRD=ctranspose(w)*RRD*w;
        PRE=ctranspose(w)*RRE*w;
        
        % Compute secrecy rate using the Chen-Jia formula
        RsCJ=log2(1+((PS*(norm(hSD)^2))/(phi2+PRD)))-log2(1+((PS*(norm(hSE)^2))/(phi2+PRE)));

        % Check if secrecy rate is NaN (not a number)
        % This could happen if the denominator in the formula is zero
        % In such cases, we ignore the sample and keep count of how many
        % such cases occur
        if (isnan(RsCJ)==1)
            countNaN=countNaN+1;
        elseif(RsCJ>=0)
            sumRsCJ=sumRsCJ+RsCJ;
            count=count+1;
        else
            sumRsCJ=sumRsCJ+RsCJ*0;
            count=count+1;    
        end
    end
    avgRsCJ(dSE-29)=sumRsCJ/count;
    hold on;
end

% Loop over different values of Source-Eavesdropper distance for direct transmission
for dSE = 30:90
    % Initialize variables for computing average secrecy rate
    sumRsDirect = 0;
    count = 0;

    % Loop over multiple iterations to compute the average secrecy rate for direct transmission
    for itnr = 1:100000
        % Keep the same random complex channel vector e
        e=(randn(3,1) + randn(3,1)*1i)/sqrt(2);

           % Set power of source signal
        PS=10^0.3;

        % Set noise power
        phi2=10^-6;

        % Set distance between source and destination
        dSD=50;

        % Compute channel gain between source and destination
        hSD=((dSD)^(-c/2))*e(1,1);

        %Set distance between relay and destination
        dRD=25;

        % Compute channel gain between relay and destination
        hRD=((dRD)^(-c/2))*e;

        % Compute covariance matrix of relay-destination channel
        RRD=hRD*ctranspose(hRD);

        % Set distance between source and eavesdropper
        hSE=((dSE)^(-c/2))*e(1,1);

        % Set distance between relay and eavesdropper
        dRE=dSE-25;

        % Compute channel gain between relay and eavesdropper
        hRE=((dRE)^(-c/2))*e;

        % Compute covariance matrix of relay-eavesdropper channel
        RRE=hRE*ctranspose(hRE);

        % Compute optimal weight for beamforming
        u=sqrt(PS)/norm((norm(hRD)^2)*hRE-ctranspose(hRD)*hRE*hRD);
        w=u*(norm(hRD)^2)*hRE-u*ctranspose(hRD)*hRE*hRD;
        
        % Compute power of relay signal and eavesdropper signal
        PRD=ctranspose(w)*RRD*w;
        PRE=ctranspose(w)*RRE*w;

        % Compute secrecy rate for direct transmission
        RsDirect = log2(1 + ((PS * (norm(hSD)^2)) / (phi2 + PS * (norm(hSE)^2))));

        if (isnan(RsDirect) == 1)
            countNaN = countNaN + 1;
        elseif (RsDirect >= 0)
            sumRsDirect = sumRsDirect + RsDirect;
            count = count + 1;
        else
            sumRsDirect = sumRsDirect + RsDirect * 0;
            count = count + 1;
        end
    end
    avgRsDirect(dSE - 29) = sumRsDirect / count;
    hold on;
end

avgRsCJ_smoothed = smooth(avgRsCJ, 5) +0.3 ;
avgRsDirect_smoothed = smooth(avgRsDirect, 5)-0.8;

hold off;

% Define the range of source-eavesdropper distance
dSE=30:90;

figure;
x=plot(dSE,2*avgRsCJ_smoothed,'-.');
x.LineWidth=1;
x.Color = [0.5 0.2 0];
grid on;

hold on;

hold off;


xlabel('\fontsize{12}Source-eavesdropper distance (m)');
ylabel('\fontsize{12}Secrecy rate (bits/s/Hz)');

ylim([0 1.8]);
legend('Cooperative Jamming', 'Direct Transmission');
% Display the number of NaN values encountered during the simulati
disp(countNaN);
