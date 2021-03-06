clc; close all; clear 
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesFontSize',14);


f = 0.2; % Porisity
E0 = 2230.6;
plotTF = false;
meshtype = {'CSTiso','LSTiso','Q4iso'};

n_max = 10;
for i = 1:length(meshtype)
    mesh = meshtype{i};
    for j = 1:n_max
        % A LST-element always has twice as many dofs as a CST or Q4. Consequently the loop
        % must run over half the nubmer of total elements if one wants to compare
        % element proformance for the same number of dofs.
        if strcmp(mesh,'LSTiso')
            nc = 2*j;
            no = nc;
        else
            nc = 4*j;
            no = nc;
        end
        
        [~,E,nu,SigmaVM_max,n_elems,eta] = FEM_master(f,mesh,nc,no,plotTF,mesh);
        E_n(j,i) = E;
        nu_n(j,i) = nu;
        SigmaVM_max_n(j,i) = SigmaVM_max;
        eta_n(j,i) = eta;       
    end
end
nc = 1:n_max;
no = nc;
dofs = 2*(2*no+1).*(2*nc+1);


%%
figure
hold on
yyaxis left
plot(dofs,E_n(:,1),'r-o')
plot(dofs,E_n(:,2),'b-^')
plot(dofs,E_n(:,3),'g-x')
xlabel('Degrees of Freedom')
ylabel('$\bar{E}$')


set(gca,{'ycolor'},{[0,0,0]})

yyaxis right
hold on
plot(dofs,nu_n(:,1),'k--o')
plot(dofs,nu_n(:,2),'m--^')
plot(dofs,nu_n(:,3),'y--x')
ylabel('Poisson ratio $\nu$')
l = legend('$\bar{E}$ CST','$\bar{E}$ LST','$\bar{E}$ Q4', '$\nu$ CST','$\nu$ LST','$\nu$ Q4','Location','East');

set(gca,{'ycolor'},{[0,0,0]})

figure
hold on
plot(dofs,SigmaVM_max_n(:,1),'r--o')
plot(dofs,SigmaVM_max_n(:,2),'b--^')
plot(dofs,SigmaVM_max_n(:,3),'g--x')
xlabel('Degrees of Freedom')
ylabel('$\sigma_{von}$')
legend('CST','LST','Q4');

figure
hold on
plot(dofs,eta_n(:,1),'r--o')
plot(dofs,eta_n(:,2),'b--^')
plot(dofs,eta_n(:,3),'g--x')
xlabel('Degrees of Freedom')
ylabel('$\eta$')
legend('CST','LST','Q4');


