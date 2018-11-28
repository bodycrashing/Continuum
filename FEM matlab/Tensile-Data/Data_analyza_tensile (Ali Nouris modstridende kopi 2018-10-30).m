clc; clear all; close all;

%% zero mm hole element
load zero_mm_one.txt
S_zero_one = zero_mm_one(:,2); %strain [mm]
F_zero_one = zero_mm_one(:,4); %force  [N] 

load zero_mm_two.txt
S_zero_two = zero_mm_two(:,2); %strain  [mm]
F_zero_two = zero_mm_two(:,4); %force   [N]

load zero_mm_three.txt
S_zero_three = zero_mm_three(:,2); %strain  [mm]
F_zero_three = zero_mm_three(:,4); %force    [N]

%Average data for zero mm holes elements
S_zero_ava = (S_zero_one+S_zero_two+S_zero_three)/3; %average strain [mm]
F_zero_ava = (F_zero_one+F_zero_two+F_zero_three)/3; %average force   [N]

%% five mm hole element
load five_mm_one.txt
S_five_one = five_mm_one(:,2); %strain [mm]
F_five_one = five_mm_one(:,4); %force   [N]

load five_mm_two.txt
S_five_two = five_mm_two(:,2); %strain  [mm]
F_five_two = five_mm_two(:,4); %force    [N]

load five_mm_three.txt
S_five_three = five_mm_three(:,2); %strain  [mm]
F_five_three = five_mm_three(:,4); %force    [N]

%Average data for five mm holes elements
S_five_ava = (S_five_one+S_five_two+S_five_three)/3; %average strain  [mm]
F_five_ava = (F_five_one+F_five_two+F_five_three)/3; %average force    [N]

%% six mm hole element
load six_mm_one.txt
S_six_one = six_mm_one(:,2); %strain  [mm]
F_six_one = six_mm_one(:,4); %force    [N]

load six_mm_two.txt
S_six_two = six_mm_two(:,2); %strain  [mm]
F_six_two = six_mm_two(:,4); %force    [N]

load six_mm_three.txt
S_six_three = six_mm_three(:,2); %strain  [mm]
F_six_three = six_mm_three(:,4); %force    [N]

%Average data for six mm holes elements
S_six_ava = (S_six_one+S_six_two+S_six_three)/3; %average strain  [mm]
F_six_ava = (F_six_one+F_six_two+F_six_three)/3; %average force    [N]

%% seven mm hole element
load seven_mm_one.txt
S_seven_one = seven_mm_one(:,2); %strain  [mm]
F_seven_one = seven_mm_one(:,4); %force    [N]
 
load seven_mm_two.txt
S_seven_two = seven_mm_two(:,2); %strain  [mm]
F_seven_two = seven_mm_two(:,4); %force    [N]

load seven_mm_three.txt
S_seven_three = seven_mm_three(:,2); %strain [mm]
F_seven_three = seven_mm_three(:,4); %force   [N]

%Average data for seven mm holes elements
S_seven_ava = (S_seven_one+S_seven_two+S_seven_three)/3; %average strain  [mm]
F_seven_ava = (F_seven_one+F_seven_two+F_seven_three)/3; %average force    [N]
%% eight mm hole element
load eight_mm_one.txt
S_eight_one = eight_mm_one(:,2); %strain  [mm]
F_eight_one = eight_mm_one(:,4); %force    [N]

load eight_mm_two.txt
S_eight_two = eight_mm_two(:,2); %strain  [mm]
F_eight_two = eight_mm_two(:,4); %force    [N]

load eight_mm_three.txt
S_eight_three = eight_mm_three(:,2); %strain  [mm]
F_eight_three = eight_mm_three(:,4); %force    [N]

%Average data for eight mm holes elements
S_eight_ava = (S_eight_one+S_eight_two+S_eight_three)/3; %average strain  [mm]
F_eight_ava = (F_eight_one+F_eight_two+F_eight_three)/3; %average force    [N]

%% plotting 
subplot(2,3,1)
plot(S_zero_one,F_zero_one);
legend('zero.mm.one tensile test')
hold on
plot(S_zero_two,F_zero_two);
legend('zero.mm.two tensile test')
hold on
plot(S_zero_three,F_zero_three);
legend('zero.mm.three tensile test')
hold on
plot(S_zero_ava,F_zero_ava);
xlim([0 1.01])
xlabel('strain')
ylabel('Force')
legend('location','northwest','zero.mm.one tensile test','zero.mm.two tensile test','zero.mm.three tensile test','zero.mm.ava tensile test')
subplot(2,3,2)
plot(S_five_one,F_five_one);
legend('five.mm.one tensile test')
hold on
plot(S_five_two,F_five_two);
legend('five.mm.two tensile test')
hold on
plot(S_five_three,F_five_three);
legend('five.mm.three tensile test')
hold on
plot(S_five_ava,F_five_ava);
xlim([0 1.01])
xlabel('strain')
ylabel('Force')
legend('location','northwest','five.mm.one tensile test','five.mm.two tensile test','five.mm.three tensile test','five.mm.ava tensile test')
subplot(2,3,3)
plot(S_six_one,F_six_one);
legend('six.mm.one tensile test')
hold on
plot(S_six_two,F_six_two);
legend('six.mm.two tensile test')
hold on
plot(S_six_three,F_six_three);
legend('six.mm.three tensile test')
hold on
plot(S_six_ava,F_six_ava);
xlim([0 1.01])
xlabel('strain')
ylabel('Force')
legend('location','northwest','six.mm.one tensile test','six.mm.two tensile test','six.mm.three tensile test','six.mm.ava tensile test')
subplot(2,3,4)
plot(S_seven_one,F_seven_one);
legend('seven.mm.one tensile test')
hold on
plot(S_seven_two,F_seven_two);
legend('seven.mm.two tensile test')
hold on
plot(S_seven_three,F_seven_three);
legend('seven.mm.three tensile test')
hold on
plot(S_seven_ava,F_seven_ava);
xlim([0 1.01])
xlabel('strain')
ylabel('Force')
legend('location','northwest','seven.mm.one tensile test','seven.mm.two tensile test','seven.mm.three tensile test','seven.mm.ava tensile test')
subplot(2,3,5)
plot(S_eight_one,F_eight_one);
legend('eight.mm.one tensile test')
hold on
plot(S_eight_two,F_eight_two);
legend('eight.mm.two tensile test')
hold on
plot(S_eight_three,F_eight_three);
legend('eight.mm.three tensile test')
hold on
plot(S_eight_ava,F_eight_ava);
xlim([0 1.01])
xlabel('strain')
ylabel('Force')
legend('location', 'northwest','eight.mm.one tensile test','eight.mm.two tensile test','eight.mm.three tensile test','eight.mm.ava tensile test')
subplot(2,3,6)
plot(S_zero_ava,F_zero_two);
hold on
plot(S_five_ava,F_five_ava);
hold on
plot(S_six_ava,F_six_ava);
hold on
plot(S_seven_ava,F_seven_ava);
hold on
plot(S_eight_ava,F_eight_ava);
xlim([0 1.01])
hold on
xlabel('strain')
ylabel('Force')
legend('location', 'northwest','zero.mm.ava tensile test','five.mm.ava tensile test','six.mm.ava tensile test','seven.mm.ava tensile test','eight.mm.ava tensile test')
%% Strain
Le= 50; % [mm] Extensometer gauge length
strian_zero_one = S_zero_one./Le;

%% Stress
Crosshead_1 = 90; %[mm^2]
Crosshead_2 = 40; %[mm^2]
Stress_zero_one= F_zero_one./Crosshead_1 ; %[N/mm^2];[MPa]
Stress_zero_one= F_zero_one./Crosshead_2 ; %[N/mm^2];[MPa]
E_zero_one = Stress_zero_one./strian_zero_one; %[N/mm^2]; [GPa]
x2= 0.0161;
y2= 0.0153;
x1= 0.0073;
y1= 0.0070;
m= (x2-x1)./(y2-y1);
E= m;
figure 
plot(strian_zero_one,Stress_zero_one)
xlabel('strain')
ylabel('Stress')
grid on
legend('Strain-Stress-line-1')
title('Stress & strain curve')

