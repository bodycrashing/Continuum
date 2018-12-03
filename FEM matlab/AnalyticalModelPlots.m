
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
ylim([0 1])
subplot(2,1,2)
plot(f,nuBar)
xlabel('Porosity')
ylabel('Poisons ratio fraction')
ylim([0 1])

%% Self-consistent model

f= 0:0.01:1;
Ebar = 1.*(1-3*f);

nuBar = 1*(1-3*f)+f;

figure
title('Dilute model')
subplot(2,1,1)
plot(f,Ebar)
xlabel('Porosity')
ylabel('Youngs modulus fraction')
ylim([0 1])
xlim([0 1])
subplot(2,1,2)
plot(f,nuBar)
xlabel('Porosity')
ylabel('Poisons ratio fraction')
ylim([0 1])
xlim([0 1])
