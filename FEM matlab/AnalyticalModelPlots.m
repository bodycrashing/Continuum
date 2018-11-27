
%% Dilute model
f= 0:0.01:1;
Ebar = 1./(1+3*f);

nuBar = (1+f)./(1+3*f);

figure
title('Dilute model')
subplot(2,1,1)
plot(f,Ebar)
xlabel('Porosity')
ylabel('Youngs modulus fraction')
subplot(2,1,2)
plot(f,nuBar)
xlabel('Porosity')
ylabel('Poisons ratio fraction')