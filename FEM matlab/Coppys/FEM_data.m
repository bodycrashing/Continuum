
close all

D = [0 2 5 6 7 8];

%porous = A_hole/

spec.E = 2.2306e+03;
spec.nu = 0.3550;

spec.D = D(2);
spec.W = 10;
spec.no = 8;
spec.nc = 16;
spec.plot = true;


[E,nu]=FEM_Solver(spec)

spec.plot = false
for i = 1:length(D)
    
    spec.D = D(i);
    E(i)=FEM_Solver(spec);
    
    
end
