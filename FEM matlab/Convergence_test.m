clear; close all

%% E_modulus test
f = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35];
for i = 1:8


[res,E_unitCell(i),nu_unitCell(i),sigma_max_Q,n_elems_Q] = FEM_master(f(i),'LSTiso',16,16,false);

end
subplot(2,1,1)
plot(f,E_unitCell)
xlabel('Porosity')
ylabel('Youngs modulus')
subplot(2,1,2)
plot(f,nu_unitCell)
xlabel('Porosity')
ylabel('Poissons ratio')
%% Convergence test

f = 0.3;
meshType = 'Q4iso'
for i = 1:10

    nc = i*4;
    no = nc
     
    [res,E_unitCell,nu_unitCell,sigma_max_Q(i),n_elems_Q(i)] = FEM_master(f,meshType,nc,no,false);
end

    meshType = 'CSTiso'
for i = 1:10

    nc = i*2;
    no = nc
     
    [res,E_unitCell,nu_unitCell,sigma_max_C(i),n_elems_C(i)] = FEM_master(f,meshType,nc,no,false);
end

    meshType = 'LSTiso'
for i = 1:6

    nc = 2*i;
    no = nc
     
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

