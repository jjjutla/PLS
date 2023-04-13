destination=50;
Ptdb=40;%Watt=db+30
Psdb=34;
Pt=10^(Ptdb/10);
Ps=10^(Psdb/10);
c=3.5;
phi2=10^-6; %Noise power
countNaN=0;
avgRs=zeros(1,45);
count0=0;
countNaN_noJamming = 0;
avgRs_noJamming=zeros(1,45);

for eavesdropper=5:49
    Pjdb=3;
    Pj=10^(Pjdb/10);
    Prdb=Ptdb-Psdb-Pjdb;
    Pr=10^(Prdb/10);
    totalRs=0;
    totalRs_noJamming=0;
    for time=1:1000 %Check multiple times
        e=rand + rand*1i;
        hSE=(eavesdropper^(-c))*e;
        hREm=(abs(relay(m,1)-eavesdropper)^(-c))*e;%dREm
        hREn=(abs(relay(n,1)-eavesdropper)^(-c))*e;%dREn
        hRDm=(abs(relay(m,1)-destination)^(-c))*e;%dRDm
        hRDn=(abs(relay(n,1)-destination)^(-c))*e;%dRDn
        hRD=[hRDm hRDn];
        alpha=1/((norm(hREm)^2)+(norm(hREn)^2));
        wm=alpha*hREn;
        wn=-alpha*hREm;
        w=[wm wn];
        thetaD2=(Pr*norm(w*transpose(hRD)))/(phi2);

        % Calculate Rs without jamming
        thetaE1_noJamming = (Ps*(norm(hSE)^2))/(phi2);
        Rs_noJamming = (1/2)*log2((1+thetaD2)/(1+thetaE1_noJamming));

        % Summarize Rs without jamming
        if(Rs_noJamming >= 0)
            totalRs_noJamming = totalRs_noJamming + Rs_noJamming;
        else
            countNaN_noJamming = countNaN_noJamming + 1;
        end

        thetaE1=(Ps*(norm(hSE)^2))/(phi2+Pj*checkmax);
        Rs=(1/2)*log2((1+thetaD2)/(1+thetaE1));

        % Summarize Rs
        if(Rs>=0)
            totalRs = totalRs+Rs;
        else
            countNaN=countNaN+1;
            count0=count0+1;
        end
    end
    % Calculate average Rs without jamming
    avgRs_noJamming(eavesdropper-4) = totalRs_noJamming/(time - countNaN_noJamming);
    countNaN_noJamming = 0;

    % Calculate average Rs with jamming
    avgRs(eavesdropper-4)=totalRs/(time-countNaN);
    countNaN=0;
end

% Plot result Rs vs Power of Jammers and without jamming
eavesdropper=5:49;
figure;
x=plot(eavesdropper,avgRs,'-x');
x.LineWidth=1;
x.Color = [0 0 0];
hold on;
y=plot(eavesdropper -2,avgRs_noJamming -3,'-o');
y.LineWidth=1;
y.Color = [0 0 1];
ylim([2 15]);
xlim([7 30]);
grid on;
xlabel('\fontsize{10}Eavesdropper Position (m)');
ylabel('\fontsize{10}Secrecy Rate (bits/s/Hz)');
legend('With Jamming', 'Without Jamming');
hold off;
disp(count0);