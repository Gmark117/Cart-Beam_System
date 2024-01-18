%%Preliminary Commands
clc;
close all;
clear all;

%%Known Quantities
m_cart = 0.3759; %Kg
m_disk = 0.1396; %Kg
m_tot = m_cart + m_disk; %Kg

%%Experimental Data Acquisition
OneDOF1 = readtable('Laboratory_Data\1dof_1.txt');
OneDOF2 = readtable('Laboratory_Data\1dof_2.txt');
OneDOF3 = readtable('Laboratory_Data\1dof_3.txt');
OneDOF4 = readtable('Laboratory_Data\1dof_4.txt');
OneDOF5 = readtable('Laboratory_Data\1dof_5.txt');

time1 = OneDOF1(:,1);
s_noise1 = OneDOF1(:,3);
time1 = table2array(time1);
s_noise1 = table2array(s_noise1);

time2 = OneDOF2(:,1);
s_noise2 = OneDOF2(:,3);
time2 = table2array(time2);
s_noise2 = table2array(s_noise2);

time3 = OneDOF3(:,1);
s_noise3 = OneDOF3(:,3);
s_noise3 = table2array(s_noise3);
time3 = table2array(time3);

time4 = OneDOF4(:,1);
s_noise4 = OneDOF4(:,3);
time4 = table2array(time4);
s_noise4 = table2array(s_noise4);

time5 = OneDOF5(:,1);
s_noise5 = OneDOF5(:,3);
time5 = table2array(time5);
s_noise5 = table2array(s_noise5);

%%Signal Filtering
n1 = 1000;
n2 = 1000;
tot = length(s_noise1);
s_temp1 = smooth(s_noise1,n1/tot);
s_smooth1 = smooth(s_temp1,n2/tot);
tot = length(s_noise2);
s_temp2 = smooth(s_noise2,n1/tot);
s_smooth2 = smooth(s_temp2,n2/tot);
tot = length(s_noise3);
s_temp3 = smooth(s_noise3,n1/tot);
s_smooth3 = smooth(s_temp3,n2/tot);
tot = length(s_noise4);
s_temp4 = smooth(s_noise4,n1/tot);
s_smooth4 = smooth(s_temp4,n2/tot);
tot = length(s_noise5);
s_temp5 = smooth(s_noise5,n1/tot);
s_smooth5 = smooth(s_temp5,n2/tot);

%%Peak Finding
[p_val1,p_loc1] = findpeaks(s_smooth1,'MinPeakDistance',0.09,'MinPeakHeight',0.1);
[p_val2,p_loc2] = findpeaks(s_smooth2,'MinPeakDistance',0.09,'MinPeakHeight',0.1);
[p_val3,p_loc3] = findpeaks(s_smooth3,'MinPeakDistance',0.09,'MinPeakHeight',0.1);
[p_val4,p_loc4] = findpeaks(s_smooth4,'MinPeakDistance',0.09,'MinPeakHeight',0.1);
[p_val5,p_loc5] = findpeaks(s_smooth5,'MinPeakDistance',0.09,'MinPeakHeight',0.1);

%%Plot First Test
figure
plot(time1,s_noise1,'Color',[0,0,1,0.1])

hold on

plot(time1,s_smooth1,'r')

hold on

plot(time1(p_loc1),s_noise1(p_loc1),'rx')
xlabel('Time')
ylabel('m/s^2')
xlim([2.5,4.5]);
ylim([-25,25]);
title('Cart Acceleration 1st Test')
legend('Unfiltered','Filtered','Peaks')

hold off
saveas(gcf, 'Plots\1. Cart Acceleration 1st Test.png');

%%Plot Second Test
figure
plot(time2,s_noise2,'Color',[0,0,1,0.1])

hold on

plot(time2,s_smooth2,'r')

hold on

plot(time2(p_loc2),s_temp2(p_loc2),'rx')
xlabel('Time')
ylabel('m/s^2')
xlim([1.4,3.4]);
ylim([-25,25]);
title('Cart Acceleration 2nd Test')
legend('Unfiltered','Filtered','Peaks')

hold off
saveas(gcf, 'Plots\2. Cart Acceleration 2nd Test.png');

%%Plot Third Test
figure
plot(time3,s_noise3,'Color',[0,0,1,0.1])

hold on

plot(time3,s_smooth3,'r')

hold on

plot(time3(p_loc3),s_temp3(p_loc3),'rx')
xlabel('Time')
ylabel('m/s^2')
xlim([1.8,3.5]);
ylim([-25,25]);
title('Cart Acceleration 3rd Test')
legend('Unfiltered','Filtered','Peaks')

hold off
saveas(gcf, 'Plots\3. Cart Acceleration 3rd Test.png');

%%Plot Fourth Test
figure
plot(time4,s_noise4,'Color',[0,0,1,0.1])

hold on

plot(time4,s_smooth4,'r')

hold on

plot(time4(p_loc4),s_temp4(p_loc4),'rx')
xlabel('Time')
ylabel('m/s^2')
xlim([1.7,3.7]);
ylim([-25,25]);
title('Cart Acceleration 4th Test')
legend('Unfiltered','Filtered','Peaks')

hold off
saveas(gcf, 'Plots\4. Cart Acceleration 4th Test.png');

%%Plot Fifth Test
figure
plot(time5,s_noise5,'Color',[0,0,1,0.1])

hold on

plot(time5,s_smooth5,'r')

hold on

plot(time5(p_loc5),s_temp5(p_loc5),'rx')
xlabel('Time')
ylabel('m/s^2')
xlim([2,3.5]);
ylim([-25,25]);
title('Cart Acceleration 5th Test')
legend('Unfiltered','Filtered','Peaks')

hold off
saveas(gcf, 'Plots\5. Cart Acceleration 5th Test.png');

%%Damping Ratio
[x_val1,x_loc1] = max(p_val1);
[x_val2,x_loc2] = max(p_val2);
[x_val3,x_loc3] = max(p_val3);
[x_val4,x_loc4] = max(p_val4);
[x_val5,x_loc5] = max(p_val5);

Delta1 = log(p_val1(x_loc1+1)/p_val1(x_loc1+2))/m_tot;
Delta2 = log(p_val2(x_loc2+1)/p_val2(x_loc2+2))/m_tot;
Delta3 = log(p_val3(x_loc3+1)/p_val3(x_loc3+2))/m_tot;
Delta4 = log(p_val4(x_loc4+1)/p_val4(x_loc4+2))/m_tot;
Delta5 = log(p_val5(x_loc5+1)/p_val5(x_loc5+2))/m_tot;

Xi1 = Delta1/sqrt(4*(pi^2)+Delta1^2);
Xi2 = Delta2/sqrt(4*(pi^2)+Delta2^2);
Xi3 = Delta3/sqrt(4*(pi^2)+Delta3^2);
Xi4 = Delta4/sqrt(4*(pi^2)+Delta4^2);
Xi5 = Delta5/sqrt(4*(pi^2)+Delta5^2);
Xi = [Xi1,Xi2,Xi3,Xi4,Xi5];

MeanXi = mean(Xi);
SigmaXi = std(Xi,1);

%%Period
T1 = mean(diff(time1(p_loc1(x_loc1+1:end)))); %s
T2 = mean(diff(time2(p_loc2(x_loc2+1:end)))); %s
T3 = mean(diff(time3(p_loc3(x_loc3+1:end)))); %s
T4 = mean(diff(time4(p_loc4(x_loc4+1:end)))); %s
T5 = mean(diff(time5(p_loc5(x_loc5+1:end)))); %s
T = [T1,T2,T3,T4,T5];

MeanT = mean(T);
SigmaT = std(T,1);

%%Damped Natural Frequency
OmegaD1 = 2*pi/T1; %Hz
OmegaD2 = 2*pi/T2; %Hz
OmegaD3 = 2*pi/T3; %Hz
OmegaD4 = 2*pi/T4; %Hz
OmegaD5 = 2*pi/T5; %Hz

%%Natural Frequency
OmegaN1 = OmegaD1/sqrt(1-(Xi1^2)); %rad/s
OmegaN2 = OmegaD2/sqrt(1-(Xi2^2)); %rad/s
OmegaN3 = OmegaD3/sqrt(1-(Xi3^2)); %rad/s
OmegaN4 = OmegaD4/sqrt(1-(Xi4^2)); %rad/s
OmegaN5 = OmegaD5/sqrt(1-(Xi5^2)); %rad/s
OmegaN = [OmegaN1,OmegaN2,OmegaN3,OmegaN4,OmegaN5];

MeanOmegaN = mean(OmegaN);
SigmaOmegaN = std(OmegaN,1);

%%Stiffness Coefficient
K1 = m_tot*(OmegaN1^2);
K2 = m_tot*(OmegaN2^2);
K3 = m_tot*(OmegaN3^2);
K4 = m_tot*(OmegaN4^2);
K5 = m_tot*(OmegaN5^2);
K = [K1,K2,K3,K4,K5];

MeanK = mean(K);
SigmaK = std(K,1);

%%Damping Coefficient
C1 = 2*m_tot*Xi1*OmegaN1;
C2 = 2*m_tot*Xi2*OmegaN2;
C3 = 2*m_tot*Xi3*OmegaN3;
C4 = 2*m_tot*Xi4*OmegaN4;
C5 = 2*m_tot*Xi5*OmegaN5;
C = [C1,C2,C3,C4,C5];

MeanC = mean(C);
SigmaC = std(C,1);

%%Data File Writing
save('Data_SingleDOF','MeanK','SigmaK','MeanC','SigmaC','MeanOmegaN','SigmaOmegaN');