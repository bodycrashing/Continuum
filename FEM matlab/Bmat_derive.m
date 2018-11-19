

%% CST

syms r s
syms x1 x2 x3 y1 y2 y3 x y

c = [x1 x2 x3;y1 y2 y3].';


N1 = 1-r-s;
N2 = r;
N3 = s;

N_CST = [N1,N2,N3];


N_diffCST = [diff(N_CST,r);diff(N_CST,s)];

J_CST = N_diffCST*c;

J_invCST = inv(J_CST);

B_CST = J_invCST*N_diffCST

%% LST
syms x1 x2 x3 x4 x5 x6 y1 y2 y3 y4 y5 y6 x y
%r = [x1 x2 x3 x4 x5 x6].';
%s = [y1 y2 y3 y4 y5 y6].';

c = [x1 x2 x3 x4 x5 x6;y1 y2 y3 y4 y5 y6].';
d = reshape(c.',1,12);

syms r s

N1 = (1-r-s)*(1-2*r-2*s);
N2 = r*(2*r-1);
N3 = s*(2*s-1);
N4 = 4*r*(1-r-s);
N5 = 4*r*s;
N6 = 4*s*(1-r-s);


N = [N1,N2,N3,N4,N5,N6];

N_diff = sym(zeros(3,12));

N_diff(1,1:2:12) = diff(N,r);
N_diff(3,2:2:12) = diff(N,r);
N_diff(2,2:2:12) = diff(N,s);
N_diff(3,1:2:12) = diff(N,s);

N_diff_square = [diff(N,r);diff(N,s)];

J = N_diff_square*c;

J_inv = inv(J);

B = J_inv*N_diff_square;

B3 = 1/det(J) * simplify(B*det(J))

simplify(B*det(J))


B2 = (1/det(J)) * N_diff;

