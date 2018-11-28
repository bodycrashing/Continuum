clc; clear all; close all;

%% zero mm hole element
load zero_mm_one.txt
S_zero_one = zero_mm_one(:,3); % Extensometer Output [mm]
F_zero_one = zero_mm_one(:,4); % Force  [N] 
%% plotting 
plot(S_zero_one,F_zero_one);
grid on
xlabel('strain')
ylabel('Force')
title('Force & strain curve')
legend('zero.mm.one tensile test')
%% Strain
Le= 50; % [mm] Extensometer gauge length
strian_zero_one = S_zero_one./Le; % Strian 
%% Stress
Crosshead_1 = 90; %[mm^2]
Stress_zero_one= F_zero_one./Crosshead_1 ; %[N/mm^2]; Stress [MPa]
figure 
plot(strian_zero_one,Stress_zero_one)
xlabel('strain')
ylabel('Stress')
grid on
legend('Strain-Stress-line-1')
title('Stress & strain curve')
x2= 0.0149;
y2= 0.0140;
x1= 0.0109;
y1= 0.0102;
m= (x2-x1)./(y2-y1);
E= m; % Elastic modulus
E_zero_one = Stress_zero_one./strian_zero_one; %[N/mm^2]; Elastic modulus [GPa]