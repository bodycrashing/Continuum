clear; close all; clc

%% E_modulus test
%f = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35];
f = 0:0.025:pi/4;
E_t  = 2230.6;
ElmType = {'LSTiso','Q4iso'};
for j = 1:length(ElmType)
    for i = 1:length(f)
    
        [res,E_unitCell(i,j),nu_unitCell(i,j),sigma_max_Q,n_elems_Q] = FEM_master(f(i),ElmType{j},18,18,false);
        
    end
end
subplot(2,1,1)
plot(f,E_unitCell(:,1)/E_t,'r',f,E_unitCell(:,2)/E_t)
xlabel('Porosity')
ylabel('Youngs modulus')
legend({ElmType{1},ElmType{2}},'FontSize',13)

subplot(2,1,2)
plot(f,nu_unitCell(:,1),'r',f,nu_unitCell(:,2))
xlabel('Porosity')
ylabel('Poissons ratio')
legend({ElmType{1},ElmType{2}},'FontSize',13)
%% Convergence test


return
f = 0.3;
meshType = 'Q4iso'
for i = 1:10

    nc = i*4;
    no = nc;
     
    [res,E_unitCell,nu_unitCell,sigma_max_Q(i),n_elems_Q(i)] = FEM_master(f,meshType,nc,no,false);
end

    meshType = 'CSTiso'
for i = 1:10

    nc = i*2;
    no = nc;
     
    [res,E_unitCell,nu_unitCell,sigma_max_C(i),n_elems_C(i)] = FEM_master(f,meshType,nc,no,false);
end

    meshType = 'LSTiso'
for i = 1:6

    nc = 2*i;
    no = nc;
     
    [res,E_unitCell,nu_unitCell,sigma_max_L(i),n_elems_L(i)] = FEM_master(f,meshType,nc,no,false);
end

%%
semilogx(1./n_elems_Q(:),sigma_max_Q(:))
hold on
semilogx(1./n_elems_C(:),sigma_max_C(:))
semilogx(1./n_elems_L(:),sigma_max_L(:))
hold off
legend('Q4','CST','LST')
xlabel('Element size')
ylabel('Max von Mises stress')

