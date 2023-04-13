clc
clear all
%%
D=2; % data relay 
N=3;  %helping relay
mc=100000;  %monte carlo simulation
eaves=3;  %number of eaves
No=(10^(-10/10))/1000;
%%
Pt_dBm=[-10:1:39];  %Global Transit power in dBm  
for j=1:length(Pt_dBm)
Ptarnsmit(j)=(10^(Pt_dBm(j)/10))/1000;  % Global Transmit power in Watt
end
%  Ptarnsmit=100;
Ptarnsmit=Ptarnsmit';
%%
%jaming signal 
jaming_signal=(10^(120/10))/1000; 
%%
%channels
G=1/sqrt(2)*((rand(D,mc)) + (rand(D,mc))*1i); %channel bewteen source and data node
h_v=1/sqrt(2)*((rand(D,mc)) + (rand(D,mc))*1i); %channel bewteen data node and vehicle
h_u=1/sqrt(2)*((rand(D,mc)) + (rand(D,mc))*1i); %channel bewteen data node and user
h_i=1/sqrt(2)*((rand(D,mc)) + (rand(D,mc))*1i); %channel bewteen data node and infrastructure
ge=1/sqrt(2)*((rand(eaves,mc)) + (rand(eaves,mc))*1i); %channel bewteen helping node and eaves
Ee=1/sqrt(2)*((rand(eaves,mc)) + (rand(eaves,mc))*1i); %channel bewteen source node and eaves
%%
theta1=angle(G);        %angle bewteen source and data node
theta2=angle(h_v);      %angle bewteen data node and vehicle
theta3=angle(h_u);     %angle bewteen data node and user
theta4=angle(h_i);      %angle bewteen data node and infrastructure
theta5=angle(ge);       %angle bewteen helping node and eaves
theta6=angle(Ee);       %angle bewteen source node and eaves
%%
% SNR at vehicle for two phase jamming 
h1=abs(h_v);
SNR_at_vehi=(h1.^2)/(No);
for i=1:mc
SNR_at_vehicle(:,i)=sum(SNR_at_vehi(:,i));
end
%%
% SNR at vehicle for simple jamming 
h1_s=(abs(h_v)).*(abs(G));
SNR_at_vehi_s=(h1_s.^2)/(No);
for i=1:mc
SNR_at_vehicle_s(:,i)=sum(SNR_at_vehi_s(:,i));
end
%%
% SNR at user for two phase jaming
h2=abs(h_u);
SNR_at_us=(h2.^2)/(No);
for i=1:mc
SNR_at_user(:,i)=sum(SNR_at_us(:,i));
end
%%
% SNR at user for jaming
h2_s=(abs(h_u)).*(abs(G));
SNR_at_us_s=(h2_s.^2)/(No);
for i=1:mc
SNR_at_user_s(:,i)=sum(SNR_at_us_s(:,i));
end
%%
% SNR at infrastructure for two phase jaming
h3=abs(h_i);
SNR_at_inf=(h3.^2)/(No);
for i=1:mc
SNR_at_infrastructure(:,i)=sum(SNR_at_inf(:,i));
end
%%
% SNR at infrastructure for jaming
h3_s=abs(h_i);
SNR_at_inf_s=(h3.^2)/(No);
for i=1:mc
SNR_at_infrastructure_s(:,i)=sum(SNR_at_inf_s(:,i));
end
%%
clear i j
% SNR at eavesdropper for two phase jamming
h4=((Ee));
h5=((ge));
theta=exp(j*theta1);
for i=1:mc
phi(:,i)=sum(theta(:,i));
end
for i=1:mc
h_4(:,i)=sum(h4(:,i));
end
for i=1:mc
h_5(:,i)=sum(h5(:,i));
end
j1=D*((((abs(h4))).*(jaming_signal)).^2); %jamming signal from data node
j2=N*(((abs(h5))*(jaming_signal)).^2); %jaming signal from helper node
channel=(abs(((h4)).*(phi))).^2;
for i=1:mc
j_1(:,i)=sum(j1(:,i));
j_2(:,i)=sum(j2(:,i));
end
for i=1:mc
Channel(:,i)=sum(channel(:,i));
end
jam=mean(j_1+j_2);
jam1=mean(j_2);
N_jam=(No)+jam;
SNR_at_eavesdropper=Channel/N_jam;
%%
clear i j
% SNR at eavesdropper for simple jaming
h4_s=((Ee));
h5_s=((ge));
theta_s=-(theta1+theta2);
for i=1:mc
phi_s(:,i)=sum(theta_s(:,i));
end
for i=1:mc
h_4_s(:,i)=sum(h4_s(:,i));
end
for i=1:mc
h_5_s(:,i)=sum(h5_s(:,i));
end
j1_s=(((abs(h5_s))*(jaming_signal)).^2); %jaming signal from helper node
channel_s=(abs(((h4_s)).*(phi_s))).^2;
for i=1:mc
j_1_s(:,i)=sum(j1_s(:,i));
end
for i=1:mc
Channel_s(:,i)=sum(channel_s(:,i));
end
jam_s=j_1_s;
N_jam_s=mean(jam_s);
SNR_at_eavesdropper_s=Channel_s/N_jam_s;
%%
Csvehicle=log2((1+SNR_at_vehicle)./(1+SNR_at_eavesdropper));
Csvehicle_s=log2((1+SNR_at_vehicle_s)./(1+SNR_at_eavesdropper_s));
Csuser=log2((1+SNR_at_user)./(1+SNR_at_eavesdropper));
Csuser_s=log2((1+SNR_at_user_s)./(1+SNR_at_eavesdropper_s));
Csinfrastructure=log2((1+SNR_at_infrastructure)./(1+SNR_at_eavesdropper));
Csinfrastructure_s=log2((1+SNR_at_infrastructure_s)./(1+SNR_at_eavesdropper_s));
for i=1:mc
Cs_vehicle(:,i)=sqrt(Ptarnsmit).*(Csvehicle(:,i));
Cs_vehicle_s(:,i)=sqrt(Ptarnsmit).*(Cs_vehicle(:,i));
Cs_user(:,i)=sqrt(Ptarnsmit).*(Csuser(:,i));
Cs_user_s(:,i)=sqrt(Ptarnsmit).*(Csuser_s(:,i));
Cs_infrastructure(:,i)=sqrt(Ptarnsmit).*(Csinfrastructure(:,i));
Cs_infrastructure_s(:,i)=sqrt(Ptarnsmit).*(Csinfrastructure_s(:,i));
end
for j=1:length(Ptarnsmit)
Rate_vehicle(j,:)=mean(Cs_vehicle(j,:));
Rate_vehicle_s(j,:)=mean(Cs_vehicle_s(j,:));
Rate_user(j,:)=mean(Cs_user(j,:));
Rate_user_s(j,:)=mean(Cs_user_s(j,:));
Rate_infrastructure(j,:)=mean(Cs_infrastructure(j,:));
Rate_infrastructure_s(j,:)=mean(Cs_infrastructure_s(j,:));
end
%%
p1 = zeros(1,length(Ptarnsmit));
p2 = zeros(1,length(Ptarnsmit));
for u=1:length(Ptarnsmit)
    for c=1:mc
      if Cs_vehicle(u,c)< 2.35
          p1(u)=p1(u)+1;
      end
    end
end
for u=1:length(Ptarnsmit)
    for c=1:mc
      if Cs_vehicle_s(u,c)< 2.35
          p2(u)=p2(u)+1;
      end
    end
end
pout=(p1/mc);
pout_c=(p2/mc);

%%
plot(Pt_dBm,pout)
hold on
plot(Pt_dBm,pout_c)
hold off
ylabel({'Secrecy Outage Probability'})
xlabel({'Transmit power'})
%plot(Pt_dBm,Rate_vehicle)
%hold on
%plot(Pt_dBm,Rate_vehicle_s)
%hold on
%plot(Pt_dBm,Rate_user,'--')
%hold on
%plot(Pt_dBm,Rate_user_s)
%hold off
%ylabel({'Secrecy Rate'})
%xlabel({'Transmit power'})
legend('Secrecy Rate of Vehicle for Two Phase Cooperative Jamming','Secrecy Rate of Vehicle for Cooperative Jamming','Secrecy Rate of User for Two Phase Cooperative Jamming','Secrecy Rate of User for Cooperative Jamming')
%legend('Secrecy rate of vehicle for two phase coperative jamming','Secrecy rate of vehicle for coperative jamming')
%legend('SOP of vehicle for two phase coperative jamming','SOP of vehicle for coperative jamming')