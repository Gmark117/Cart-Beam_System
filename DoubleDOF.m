%%Preliminary Commands
clc;
close all;
clear all;

%%Known Quantities
Young_modulus = 210; %GPa
Density = 7850; %Kg/m^3

m_cart = 0.3759; %Kg
m_disk = 0.1396; %Kg
m_beam = 4.7764; %Kg
m_shaker = 0.2000; %Kg
m_tot = m_cart + m_disk; %Kg

rod_length = 100.0; %mm
rod_width = 25.0; %mm
rod_thickness = 1.5; %mm
rod_damping_ratio = 0.01;
rod_inertia = rod_width*(rod_thickness^3)/12; %mm^4
rod_stiffness = 1e6*12*Young_modulus*rod_inertia/(rod_length^3);

m_rod = Density*rod_width*rod_thickness*rod_length*1e-9; %Kg
rod_damping_factor = 2*rod_damping_ratio*sqrt(rod_stiffness*m_rod);

beam_length = 0.605; %mm
beam_width = 0.030; %mm
beam_thickness = 0.030; %mm

Fs = 6400; %Hz
min_Freq = 5; %Hz
max_Freq = 30; %Hz
Duration = 40; %s

load('Data_SingleDOF');

%%Experimental Data Acquisition
TwoDOF1 = readtable('Laboratory_Data\2dof_1.txt');
TwoDOF2 = readtable('Laboratory_Data\2dof_2.txt');
TwoDOF3 = readtable('Laboratory_Data\2dof_3.txt');
TwoDOF4 = readtable('Laboratory_Data\2dof_4.txt');
TwoDOF5 = readtable('Laboratory_Data\2dof_5.txt');

time1 = TwoDOF1(:,1);
force1 = TwoDOF1(:,2);
cart_noise1 = TwoDOF1(:,3);
beam_noise1 = TwoDOF1(:,4);
time1 = table2array(time1);
force1 = table2array(force1);
cart_noise1 = table2array(cart_noise1);
beam_noise1 = table2array(beam_noise1);

time2 = TwoDOF2(:,1);
force2 = TwoDOF2(:,2);
cart_noise2 = TwoDOF2(:,3);
beam_noise2 = TwoDOF2(:,4);
time2 = table2array(time2);
force2 = table2array(force2);
cart_noise2 = table2array(cart_noise2);
beam_noise2 = table2array(beam_noise2);

time3 = TwoDOF3(:,1);
force3 = TwoDOF3(:,2);
cart_noise3 = TwoDOF3(:,3);
beam_noise3 = TwoDOF3(:,4);
time3 = table2array(time3);
force3 = table2array(force3);
cart_noise3 = table2array(cart_noise3);
beam_noise3 = table2array(beam_noise3);

time4 = TwoDOF4(:,1);
force4 = TwoDOF4(:,2);
cart_noise4 = TwoDOF4(:,3);
beam_noise4 = TwoDOF4(:,4);
time4 = table2array(time4);
force4 = table2array(force4);
cart_noise4 = table2array(cart_noise4);
beam_noise4 = table2array(beam_noise4);

time5 = TwoDOF5(:,1);
force5 = TwoDOF5(:,2);
cart_noise5 = TwoDOF5(:,3);
beam_noise5 = TwoDOF5(:,4);
time5 = table2array(time5);
force5 = table2array(force5);
cart_noise5 = table2array(cart_noise5);
beam_noise5 = table2array(beam_noise5);

%%POINT 1-2
h11 = tf([m_tot MeanC MeanK 0 0], ...
          [(m_beam*m_tot) ...
          (MeanC*m_beam + MeanC*m_tot + 2*rod_damping_factor*m_tot) ...
          (m_beam*MeanK + 2*rod_damping_factor*MeanC + 2*rod_stiffness*m_tot + MeanK*m_tot) ...
          (2*rod_damping_factor*MeanK + 2*rod_stiffness*MeanC) ...
          (2*rod_stiffness*MeanK)]);

h21 = tf([0 MeanC MeanK 0 0], ...
          [(m_beam*m_tot) ...
          (MeanC*m_beam + MeanC*m_tot + 2*rod_damping_factor*m_tot) ...
          (m_beam*MeanK + 2*rod_damping_factor*MeanC + 2*rod_stiffness*m_tot + MeanK*m_tot) ...
          (2*rod_damping_factor*MeanK + 2*rod_stiffness*MeanC) ...
          (2*rod_stiffness*MeanK)]);
%{
h12 = tf([0,MeanC,MeanK,0,0], ...
          [m_beam*m_tot ...
          (MeanC*m_beam + MeanC*m_tot + 2*rod_damping_factor*m_tot) ...
          (m_beam*MeanK + 2*rod_damping_factor*MeanC + 2*rod_stiffness*m_tot + MeanK*m_tot) ...
          (2*rod_damping_factor*MeanK + 2*rod_stiffness*MeanC) ...
          2*rod_stiffness*MeanK]);

h22 = tf([m_beam,MeanC + 2*rod_damping_factor,MeanK + 2*rod_stiffness,0,0], ...
          [m_beam*m_tot ...
          (MeanC*m_beam + MeanC*m_tot + 2*rod_damping_factor*m_tot) ...
          (m_beam*MeanK + 2*rod_damping_factor*MeanC + 2*rod_stiffness*m_tot + MeanK*m_tot) ...
          (2*rod_damping_factor*MeanK + 2*rod_stiffness*MeanC) ...
          2*rod_stiffness*MeanK]);
%}

%====================Error Propagation========================
%h11
errC = 3*SigmaC;
errK = 3*SigmaK;
Num_h11 = m_tot+MeanC+MeanK;
Den_h11 = m_beam*m_tot + MeanC*m_beam + MeanC*m_tot + 2*rod_damping_factor*m_tot ...
      +m_beam*MeanK + 2*rod_damping_factor*MeanC + 2*rod_stiffness*m_tot ...
      +MeanK*m_tot + 2*rod_damping_factor*MeanK + 2*rod_stiffness*MeanC ...
      +2*rod_stiffness*MeanK;
eNh11 = sqrt(errC^2+errK^2);
eDh11 = sqrt(((m_beam+m_tot+2*rod_damping_factor+2*rod_stiffness)^2)*(errC^2) ...
               +((m_beam+m_tot+2*rod_damping_factor+2*rod_stiffness)^2)*(errK^2));
eh11 = (Num_h11/Den_h11)*sqrt((eNh11/(Num_h11))^2+(eDh11/(Den_h11))^2);

h11_Upper = h11+eh11;
h11_Lower = h11-eh11;

%h21
errC = 3*SigmaC;
errK = 3*SigmaK;
Num_h21 = MeanC+MeanK;
Den_h21 = m_beam*m_tot + MeanC*m_beam + MeanC*m_tot + 2*rod_damping_factor*m_tot ...
          +m_beam*MeanK + 2*rod_damping_factor*MeanC + 2*rod_stiffness*m_tot + MeanK*m_tot ...
          +2*rod_damping_factor*MeanK + 2*rod_stiffness*MeanC ...
          +2*rod_stiffness*MeanK;
eNh21 = sqrt(errC^2+errK^2);
eDh21 = sqrt(((m_beam+m_tot+2*rod_damping_factor+2*rod_stiffness)^2)*(errC^2) ...
               +((m_beam+m_tot+2*rod_damping_factor+2*rod_stiffness)^2)*(errK^2));
eh21 = (Num_h21/Den_h21)*sqrt((eNh21/(Num_h21))^2+(eDh21/(Den_h21))^2);

h21_Upper = h21+eh21;
h21_Lower = h21-eh21;

%========================Plots================================
opts = bodeoptions('cstprefs');
opts.PhaseVisible = 'off';
opts.MagUnits = 'abs';

figure(1)
hold on
bodeplot(h11,'k',h11_Upper,'r:',h11_Lower,'r--',{4*2*pi,31*2*pi},opts);
bodeplot(h21,'k',h21_Upper,'b:',h21_Lower,'b--',{4*2*pi,31*2*pi},opts);

grid on
title('Cart/Beam Acceleration vs Force Transfer Function Bode Plot')
legend('Beam Acc. vs Force','Beam Acc. vs Force Upper','Beam Acc. vs Force Lower', ...
       'Cart Acc. vs Force','Cart Acc. vs Force Upper','Cart Acc. vs Force Lower')
%hold off
saveas(gcf, 'Plots\6. Cart-Beam Acceleration vs Force Transfer Function Bode Plot.png');

%%POINT 3
%=======================Cart==================================
[Tf1c,Fr1c] = tfestimate(force1,cart_noise1,[],[],[],Fs);
[Tf2c,Fr2c] = tfestimate(force2,cart_noise2,[],[],[],Fs);
[Tf3c,Fr3c] = tfestimate(force3,cart_noise3,[],[],[],Fs);
[Tf4c,Fr4c] = tfestimate(force4,cart_noise4,[],[],[],Fs);
[Tf5c,Fr5c] = tfestimate(force5,cart_noise5,[],[],[],Fs);
Tfc = [Tf1c,Tf2c,Tf3c,Tf4c,Tf5c];

Tm1c = abs(Tf1c);
Tm2c = abs(Tf2c);
Tm3c = abs(Tf3c);
Tm4c = abs(Tf4c);
Tm5c = abs(Tf5c);

lim1 = find((interp1(4:31,Fr1c,'nearest') == 4),1,'first');
lim2 = find((interp1(4:31,Fr1c,'nearest') == 31),1,'last');

sys1c = frd(Tf1c(lim1:lim2),Fr1c(lim1:lim2).*(2*pi));
sys2c = frd(Tf2c(lim1:lim2),Fr2c(lim1:lim2).*(2*pi));
sys3c = frd(Tf3c(lim1:lim2),Fr3c(lim1:lim2).*(2*pi));
sys4c = frd(Tf4c(lim1:lim2),Fr4c(lim1:lim2).*(2*pi));
sys5c = frd(Tf5c(lim1:lim2),Fr5c(lim1:lim2).*(2*pi));

%=========================Beam================================
[Tf1b,Fr1b] = tfestimate(force1,beam_noise1,[],[],[],Fs);
[Tf2b,Fr2b] = tfestimate(force2,beam_noise2,[],[],[],Fs);
[Tf3b,Fr3b] = tfestimate(force3,beam_noise3,[],[],[],Fs);
[Tf4b,Fr4b] = tfestimate(force4,beam_noise4,[],[],[],Fs);
[Tf5b,Fr5b] = tfestimate(force5,beam_noise5,[],[],[],Fs);
Tfb = [Tf1b,Tf2b,Tf3b,Tf4b,Tf5b];

Tm1b = abs(Tf1b);
Tm2b = abs(Tf2b);
Tm3b = abs(Tf3b);
Tm4b = abs(Tf4b);
Tm5b = abs(Tf5b);

lim3 = find((interp1(4:31,Fr1b,'nearest') == 4),1,'first');
lim4 = find((interp1(4:31,Fr1b,'nearest') == 31),1,'last');

sys1b = frd(Tf1b(lim3:lim4),Fr1b(lim3:lim4).*(2*pi));
sys2b = frd(Tf2b(lim3:lim4),Fr2b(lim3:lim4).*(2*pi));
sys3b = frd(Tf3b(lim3:lim4),Fr3b(lim3:lim4).*(2*pi));
sys4b = frd(Tf4b(lim3:lim4),Fr4b(lim3:lim4).*(2*pi));
sys5b = frd(Tf5b(lim3:lim4),Fr5b(lim3:lim4).*(2*pi));

%========================Plots================================
figure(2)
hold on
bode(sys1b,{31,190},opts);
bode(sys1c,{31,190},opts);
grid on
title('1st Test Bode Plot');
legend('Beam Acc. vs Force','Cart Acc. vs Force');
hold off
saveas(gcf, 'Plots\7. 1st Test Bode Plot.png');

figure(3)
hold on
bode(sys2b,{31,190},opts);
bode(sys2c,{31,190},opts);
grid on
title('2nd Test Bode Plot');
legend('Beam Acc. vs Force','Cart Acc. vs Force');
hold off
saveas(gcf, 'Plots\8. 2nd Test Bode Plot.png');

figure(4)
hold on
bode(sys3b,{31,190},opts);
bode(sys3c,{31,190},opts);
grid on
title('3rd Test Bode Plot');
legend('Beam Acc. vs Force','Cart Acc. vs Force');
hold off
saveas(gcf, 'Plots\9. 3rd Test Bode Plot.png');

figure(5)
hold on
bode(sys4b,{31,190},opts);
bode(sys4c,{31,190},opts);
grid on
title('4th Test Bode Plot');
legend('Beam Acc. vs Force','Cart Acc. vs Force');
hold off
saveas(gcf, 'Plots\10. 4th Test Bode Plot.png');

figure(6)
hold on
bode(sys5b,{31,190},opts);
bode(sys5c,{31,190},opts);
grid on
title('5th Test Bode Plot');
legend('Beam Acc. vs Force','Cart Acc. vs Force');
hold off
saveas(gcf, 'Plots\11. 5th Test Bode Plot.png');

%%POINT 4
%=======================Cart==================================
MeanTfc = mean(Tfc,2);
SigmaTfc_R = std(real(Tfc),0,2);
SigmaTfc_I = std(imag(Tfc),0,2);
SigmaTfc = complex(SigmaTfc_R,SigmaTfc_I);
Upper_Tfc = MeanTfc+3*SigmaTfc;
Lower_Tfc = MeanTfc-3*SigmaTfc;

sys_meanC = frd(MeanTfc(lim1:lim2),Fr1c(lim1:lim2).*(2*pi));
sys_upperC = frd(Upper_Tfc(lim1:lim2),Fr1c(lim1:lim2).*(2*pi));
sys_lowerC = frd(Lower_Tfc(lim1:lim2),Fr1c(lim1:lim2).*(2*pi));

%=======================Beam==================================
MeanTfb = mean(Tfb,2);
SigmaTfb_R = std(real(Tfb),0,2);
SigmaTfb_I = std(imag(Tfb),0,2);
SigmaTfb = complex(SigmaTfb_R,SigmaTfb_I);
Upper_Tfb = MeanTfb+3*SigmaTfb;
Lower_Tfb = MeanTfb-3*SigmaTfb;

sys_meanB = frd(MeanTfb(lim3:lim4),Fr1b(lim3:lim4).*(2*pi));
sys_upperB = frd(Upper_Tfb(lim3:lim4),Fr1b(lim3:lim4).*(2*pi));
sys_lowerB = frd(Lower_Tfb(lim3:lim4),Fr1b(lim3:lim4).*(2*pi));

%========================Plots================================
figure(7)
hold on
bodeplot(sys_meanB,'r',sys_upperB,'r:',sys_lowerB,'r--',{31,190},opts);
bodeplot(sys_meanC,'b',sys_upperC,'b:',sys_lowerC,'b--',{31,190},opts);
grid on
title('Mean Experimental Transfer Function Bode Plots');
legend('Beam Acc. vs Force', ...
       'Beam Acc. vs Force Upper', ...
       'Beam Acc. vs Force Lower', ...
       'Cart Acc. vs Force', ...
       'Cart Acc. vs Force Upper', ...
       'Cart Acc. vs Force Lower');
hold off
saveas(gcf, 'Plots\12. Mean Experimental Transfer Function Bode Plots.png');

%%POINT 5
figure(8)
hold on
bodeplot(h11,'b',h11_Upper,'b:',h11_Lower,'b--',{30,190},opts);
bodeplot(sys_meanB,'r',sys_upperB,'r:',sys_lowerB,'r--',{30,190},opts);
grid on
title('Analytical and experimental Beam Acc. vs Force Bode Plots');
legend('Analytical Beam Acc. vs Force', ...
       'Analytical Beam Acc. vs Force Upper', ...
       'Analytical Beam Acc. vs Force Lower', ...
       'Experimental Beam Acc. vs Force', ...
       'Experimental Beam Acc. vs Force Upper', ...
       'Experimental Beam Acc. vs Force Lower');
hold off
saveas(gcf, 'Plots\13. Analytical and experimental Beam Acc. vs Force Bode Plots.png');

figure(9)
hold on
bodeplot(h21,'b',h21_Upper,'b:',h21_Lower,'b--',{30,190},opts);
bodeplot(sys_meanC,'r',sys_upperC,'r:',sys_lowerC,'r--',{30,190},opts);
grid on
title('Analytical and experimental Cart Acc. vs Force Bode Plots');
legend('Analytical Cart Acc. vs Force', ...
       'Analytical Cart Acc. vs Force Upper', ...
       'Analytical Cart Acc. vs Force Lower', ...
       'Experimental Cart Acc. vs Force', ...
       'Experimental Cart Acc. vs Force Upper', ...
       'Experimental Cart Acc. vs Force Lower');
hold off
saveas(gcf, 'Plots\14. Analytical and experimental Cart Acc. vs Force Bode Plots.png');

%%POINT 6
s21 = sqrt(-1)*Fr1c(lim1:lim2).*(2*pi);
s11 = sqrt(-1)*Fr1b(lim3:lim4).*(2*pi);

par_h11 = @(K,C,rod_K,rod_C) ...
            (m_tot.*s11.^4 + C.*s11.^3 + K.*s11.^2)./ ...
            ((m_beam*m_tot).*s11.^4+ ...
            (C*m_beam + C*m_tot + 2*rod_C*m_tot).*s11.^3+ ...
            (m_beam*K + 2*rod_C*C + 2*rod_K*m_tot + K*m_tot).*s11.^2+ ...
            (2*rod_C*K + 2*rod_K*C).*s11+ ...
            2*rod_K*K);
par_h21 = @(K,C,rod_K,rod_C) ...
            (C.*s21.^3 + K.*s21.^2)./ ...
            ((m_beam*m_tot).*s21.^4+ ...
            (C*m_beam + C*m_tot + 2*rod_C*m_tot).*s21.^3+ ...
            (m_beam*K + 2*rod_C*C + 2*rod_K*m_tot + K*m_tot).*s21.^2+ ...
            (2*rod_C*K + 2*rod_K*C).*s21+ ...
            2*rod_K*K);

mod_par_h11 = abs(par_h11(2000,1,2e-05,5e-10));
mod_par_h21 = abs(par_h21(2000,1,2e-05,5e-10));

MeanTmC = abs(MeanTfc(lim1:lim2).*(2*pi));
MeanTmB = abs(MeanTfb(lim3:lim4).*(2*pi));

err_h11 = @(x) ...
            rms(MeanTmB-abs(par_h11(x(1),x(2),x(3),x(4))));
err_h21 = @(x) ...
            rms(MeanTmC-abs(par_h21(x(1),x(2),x(3),x(4))));

K0 = 2000;
C0 = 1;
rod_K0 = 1e4;
rod_C0 = 0.7;
x0 = [K0,C0,rod_K0,rod_C0];

x_opt_h11 = fminsearch(err_h11,x0);
x_opt_h21 = fminsearch(err_h21,x0);

%========================Plots================================
figure(10)
hold on
plot(Fr1c(lim1:lim2).*(2*pi),MeanTmB)
plot(Fr1c(lim1:lim2).*(2*pi),abs(par_h11(x0(1),x0(2),x0(3),x0(4))))
plot(Fr1c(lim1:lim2).*(2*pi),abs(par_h11(x_opt_h11(1),x_opt_h11(2),x_opt_h11(3),x_opt_h11(4))))
grid on
xlabel('Frequency (rad/s)')
ylabel('Module Of The Transfer Function')
title('Beam Transfer Function Fitting')
legend('Module Of Experimental Transfer Function', ...
       'Module Of Analytical Transfer Function at Initial Guess', ...
       'Module Of Analytical Transfer Function at Optimal Solution')
hold off
saveas(gcf, 'Plots\15. Beam Transfer Function Fitting.png');

figure(11)
hold on
plot(Fr1c(lim1:lim2).*(2*pi),MeanTmC)
plot(Fr1c(lim1:lim2).*(2*pi),abs(par_h21(x0(1),x0(2),x0(3),x0(4))))
plot(Fr1c(lim1:lim2).*(2*pi),abs(par_h21(x_opt_h21(1),x_opt_h21(2),x_opt_h21(3),x_opt_h21(4))))
grid on
xlabel('Frequency (rad/s)')
ylabel('Module Of The Transfer Function')
title('Cart Transfer Function Fitting')
legend('Module Of Experimental Transfer Function', ...
       'Module Of Analytical Transfer Function at Initial Guess', ...
       'Module Of Analytical Transfer Function at Optimal Solution')
hold off
saveas(gcf, 'Plots\16. Cart Transfer Function Fitting.png');










