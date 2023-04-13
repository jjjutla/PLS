n=5;
eavesdropper=25;
destination=50;
c=3.5/2;
a=1;
b=49;
relay=[3;30;33;35;37];

% Added coordinates for the specified plot
relayx = [30, 35, 40, 45, 50];
relayy = [0.059, 0.185, 0.430, 0.620, 0.700];

Rsosj=zeros(1,11);
Rsos=zeros(1,11);
Rsoscj=zeros(1,11);
totalosj=zeros(1,11);
totalos=zeros(1,11);
totaloscj=zeros(1,11);
for times=1:10000
    e=rand + rand*1i;
    for Pdb=0:5:50
        P=10^(Pdb/10);
        for m=1:n
            for o=1:n
                if o~=m
                   hRD=(abs(relay(m,1)-destination)^(-c));
                   hRE=(abs(relay(m,1)-eavesdropper)^(-c));
                   hJD=(abs(relay(o,1)-destination)^(-c));
                   hJE=(abs(relay(o,1)-eavesdropper)^(-c));
                   thetaRD=P*(norm(hRD)^2);
                   thetaRE=P*(norm(hRE)^2);
                   thetaJD=P*(norm(hJD)^2);
                   thetaJE=P*(norm(hJE)^2);

                   Rdosj=(1/2)*log2(1+(thetaRD/(1+thetaJD)));  %OSJ
                   Reosj=(1/2)*log2(1+(thetaRE/(1+thetaJE)));
                   Rsosj1=Rdosj-Reosj;
                   if Rsosj(1,Pdb/5+1)<Rsosj1
                   Rsosj(1,Pdb/5+1)=Rsosj1;
                   chosenjammerosj=o;
                   chosenrelayosj=m;
                   end

                   if thetaJE>thetaJD                           %OS
                       Rdos=(1/2)*log2(1+(thetaRD/thetaJD));  
                       Reos=(1/2)*log2(1+(thetaRE/thetaJE));
                       Rsos1=Rdos-Reos;
                       if Rsos(1,Pdb/5+1)<Rsos1
                       Rsos(1,Pdb/5+1)=Rsos1;
                       chosenjammeros=o;
                       chosenrelayos=m;
                       end 
                   end

                   Rdoscj=(1/2)*log2(1+thetaRD);
                   Reoscj=(1/2)*log2(1+(thetaRE/(1+thetaJE)));
                   Rsoscj1=Rdoscj-Reoscj;
                   if Rsoscj(1,Pdb/5+1)<Rsoscj1
                   Rsoscj(1,Pdb/5+1)=Rsoscj1;
                   chosenjammeroscj=o;
                   chosenrelayoscj=m;
                   end
                end
            end
        end    
    end
totaloscj=totaloscj+Rsoscj;
totalos=totalos+Rsos;
totalosj=totalosj+Rsosj;
end
totaloscj=totaloscj/times;
totalos=totalos/times;
totalosj=totalosj/times;
Pdb=0:5:50;

% Plot commands
x=plot(Pdb ,totalosj,'-x');
x.LineWidth=1;
hold on;
y=plot(relayx, relayy,'-o');
y.LineWidth=1;
z=plot(Pdb,totaloscj,'--');
z.LineWidth=1;
z.Color = [0 0 0];

% Plotting specified coordinates
%p = plot(specified_coordinates_x, specified_coordinates_y, 'r*');
%p.MarkerSize = 8;

% Add title, labels, and legend to the plot
grid on;
xlabel('\fontsize{10} P [dB]');
ylabel('\fontsize{10} Secrecy Ergodic Capacity');
legend('OSJ', 'OS', 'OSCJ');
hold off;

disp('Chosen relay and jammer nodes for OSJ strategy:');
disp([chosenrelayosj; chosenjammerosj]);

disp('Chosen relay and jammer nodes for OS strategy:');
disp([chosenrelayos; chosenjammeros]);

disp('Chosen relay and jammer nodes for OSCJ strategy:');
disp([chosenrelayoscj; chosenjammeroscj]);

% Display the end of the simulation
disp('end');
