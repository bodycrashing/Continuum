


f = linspace(0,pi/4,9)

f_ans = linspace(0,pi/4,5)

r = 2*1*sqrt(f/(pi))

f(end) = f(end)-0.01
s_max = [2.8e-12 5.393 4.672 4.245 4.005 3.953 4.1884 5.093 10.79]*1e+10

for i = 1:length(f)
 [res,E_unitCell,nu_unitCell,sigma_max(i)] = FEM_master(f(i),'Q4iso',36,36,false);
end

figure(1)
plot(f,s_max)
hold on
plot(f,sigma_max*1e7)
hold off
legend('Ansys','Matlab')

diff = (s_max-sigma_max*1e7)./s_max
figure(2)
plot(f,diff)

