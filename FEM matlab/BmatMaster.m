function [B,J_det,N] = BmatMaster(type,x,y,xi,eta)
% type = string containing type of element
% values of type can be: Q4iso Q8iso ....
%  xi and eta = points to evaluate B at. Must be between -1 and 1
%  x_vec and y_vec = Vertical vectors containing coordinates of element nodes
%

switch type
    case 'Q4iso'
        if nargin < 4
            error('Specify xi and eta for Q4iso element')
        end
        % Calculating B matrix for isoparametric Q4 following principles from
        % Cook p. 205-208.
        % xi and eta = points to evaluate B at. Must be between -1 and 1
        % x_vec and y_vec = Vertical vectors containing coordinates of element nodes
        
        % Calculating the Jacobian matrix, eq 6.2-6 in cook
        J1 = 1/4 * [-(1-eta) (1-eta) (1+eta) -(1+eta);
            -(1-xi) -(1+xi)  (1+xi)   (1-xi)];
        
        J = J1*cat(2,x,y);
        % Determinant of J
        J_det = det(J);
        % Finding inverse of J as this is needed to calculate B
        I = inv(J);
        
        % Following the procedure from p. 208 in cook.
        % Multiplying together eq. 6.2-9 * 6.2-10 * 6.2-11 gives us the stress
        % strain relations, which is the B matrix.
        B1 = [1 0 0 0;
            0 0 0 1;
            0 1 1 0] ;
        
        B2 = blkdiag(I,I);
        
        B3 = 1/4 * [-1+eta 0 1-eta 0 1+eta 0 -1-eta 0;
            -1+xi  0 -1-xi 0 1+xi  0  1-xi  0;
            0 -1+eta 0 1-eta 0 1+eta 0 -1-eta;
            0 -1+xi  0 -1-xi 0 1+xi  0  1-xi];
        
        B = B1*B2*B3;
        
        N = [];
    case 'Q8iso'
        
    case 'CSTiso'
        % from Cook chap 7.2
        A = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        
        B = 1/(2*A) * [y(2)-y(3), 0 ,y(3)-y(1), 0, y(1)-y(2),0;
                       0, x(3)-x(2), 0 ,x(1)-x(3), 0, x(2)-x(1);
                       x(3)-x(2), y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];
       J_det = x(1)*y(2) - x(2)*y(1) - x(1)*y(3) + x(3)*y(1) + x(2)*y(3) - x(3)*y(2);
                  N = []
    case 'LSTiso'
        r = xi;
        s = eta;
        
        N = [(2*r + 2*s - 1)*(r + s - 1),0, r*(2*r - 1), 0 , s*(2*s - 1),0, -4*r*(r + s - 1),0, 4*r*s,0, -4*s*(r + s - 1),0;
             0, (2*r + 2*s - 1)*(r + s - 1),0, r*(2*r - 1), 0 , s*(2*s - 1),0, -4*r*(r + s - 1),0, 4*r*s,0, -4*s*(r + s - 1)];
        J_det = 3*x(2)*y(1) - 3*x(1)*y(2) + 3*x(1)*y(3) - 3*x(3)*y(1) + 12*x(1)*y(4) + x(2)*y(3) - x(3)*y(2) - 12*x(4)*y(1) - 12*x(1)*y(6) + 4*x(3)*y(4) - 4*x(4)*y(3) + 12*x(6)*y(1) - 4*x(2)*y(6) + 4*x(6)*y(2) + 16*x(4)*y(6) - 16*x(6)*y(4) - 16*r^2*x(1)*y(2) + 16*r^2*x(2)*y(1) + 16*r^2*x(1)*y(4) - 16*r^2*x(4)*y(1) + 16*r^2*x(1)*y(5) - 16*r^2*x(2)*y(4) + 16*r^2*x(4)*y(2) - 16*r^2*x(5)*y(1) - 16*r^2*x(1)*y(6) + 16*r^2*x(2)*y(5) - 16*r^2*x(5)*y(2) + 16*r^2*x(6)*y(1) - 16*r^2*x(2)*y(6) + 16*r^2*x(6)*y(2) - 32*r^2*x(4)*y(5) + 32*r^2*x(5)*y(4) + 32*r^2*x(4)*y(6) - 32*r^2*x(6)*y(4) + 16*s^2*x(1)*y(3) - 16*s^2*x(3)*y(1) + 16*s^2*x(1)*y(4) - 16*s^2*x(4)*y(1) - 16*s^2*x(1)*y(5) + 16*s^2*x(5)*y(1) - 16*s^2*x(1)*y(6) + 16*s^2*x(3)*y(4) - 16*s^2*x(4)*y(3) + 16*s^2*x(6)*y(1) - 16*s^2*x(3)*y(5) + 16*s^2*x(5)*y(3) + 16*s^2*x(3)*y(6) - 16*s^2*x(6)*y(3) + 32*s^2*x(4)*y(6) - 32*s^2*x(6)*y(4) - 32*s^2*x(5)*y(6) + 32*s^2*x(6)*y(5) + 16*r*x(1)*y(2) - 16*r*x(2)*y(1) - 4*r*x(1)*y(3) + 4*r*x(3)*y(1) - 28*r*x(1)*y(4) - 4*r*x(2)*y(3) + 4*r*x(3)*y(2) + 28*r*x(4)*y(1) - 12*r*x(1)*y(5) + 4*r*x(2)*y(4) - 4*r*x(4)*y(2) + 12*r*x(5)*y(1) + 28*r*x(1)*y(6) - 4*r*x(2)*y(5) - 8*r*x(3)*y(4) + 8*r*x(4)*y(3) + 4*r*x(5)*y(2) - 28*r*x(6)*y(1) + 20*r*x(2)*y(6) - 20*r*x(6)*y(2) + 16*r*x(4)*y(5) - 16*r*x(5)*y(4) - 48*r*x(4)*y(6) + 48*r*x(6)*y(4) + 4*s*x(1)*y(2) - 4*s*x(2)*y(1) - 16*s*x(1)*y(3) + 16*s*x(3)*y(1) - 28*s*x(1)*y(4) - 4*s*x(2)*y(3) + 4*s*x(3)*y(2) + 28*s*x(4)*y(1) + 12*s*x(1)*y(5) - 12*s*x(5)*y(1) + 28*s*x(1)*y(6) - 20*s*x(3)*y(4) + 20*s*x(4)*y(3) - 28*s*x(6)*y(1) + 8*s*x(2)*y(6) + 4*s*x(3)*y(5) - 4*s*x(5)*y(3) - 8*s*x(6)*y(2) - 4*s*x(3)*y(6) + 4*s*x(6)*y(3) - 48*s*x(4)*y(6) + 48*s*x(6)*y(4) + 16*s*x(5)*y(6) - 16*s*x(6)*y(5) - 16*r*s*x(1)*y(2) + 16*r*s*x(2)*y(1) + 16*r*s*x(1)*y(3) - 16*r*s*x(3)*y(1) + 32*r*s*x(1)*y(4) + 16*r*s*x(2)*y(3) - 16*r*s*x(3)*y(2) - 32*r*s*x(4)*y(1) - 32*r*s*x(1)*y(6) + 32*r*s*x(3)*y(4) - 32*r*s*x(4)*y(3) + 32*r*s*x(6)*y(1) - 32*r*s*x(2)*y(6) + 32*r*s*x(6)*y(2) + 64*r*s*x(4)*y(6) - 64*r*s*x(6)*y(4);
        
        B_SF = 1/J_det * [  (4*r + 4*s - 3)*(y(2) - y(3) - 4*y(4) + 4*y(6) - 4*r*y(2) + 4*r*y(4) + 4*r*y(5) - 4*r*y(6) + 4*s*y(3) + 4*s*y(4) - 4*s*y(5) - 4*s*y(6)), -(4*r - 1)*(3*y(1) + y(3) - 4*y(6) - 4*r*y(1) + 4*r*y(4) - 4*r*y(5) + 4*r*y(6) - 4*s*y(1) - 4*s*y(3) + 8*s*y(6)),  (4*s - 1)*(3*y(1) + y(2) - 4*y(4) - 4*r*y(1) - 4*r*y(2) + 8*r*y(4) - 4*s*y(1) + 4*s*y(4) - 4*s*y(5) + 4*s*y(6)), 16*y(6) - 4*y(3) - 12*y(1) + 28*r*y(1) - 4*r*y(2) + 8*r*y(3) + 16*r*y(5) - 48*r*y(6) + 28*s*y(1) + 20*s*y(3) - 48*s*y(6) - 16*r^2*y(1) + 16*r^2*y(2) - 32*r^2*y(5) + 32*r^2*y(6) - 16*s^2*y(1) - 16*s^2*y(3) + 32*s^2*y(6) - 32*r*s*y(1) - 32*r*s*y(3) + 64*r*s*y(6), 12*r*y(1) + 4*r*y(2) - 16*r*y(4) - 12*s*y(1) - 4*s*y(3) + 16*s*y(6) - 16*r^2*y(1) - 16*r^2*y(2) + 32*r^2*y(4) + 16*s^2*y(1) + 16*s^2*y(3) - 32*s^2*y(6), 12*y(1) + 4*y(2) - 16*y(4) - 28*r*y(1) - 20*r*y(2) + 48*r*y(4) - 28*s*y(1) - 8*s*y(2) + 4*s*y(3) + 48*s*y(4) - 16*s*y(5) + 16*r^2*y(1) + 16*r^2*y(2) - 32*r^2*y(4) + 16*s^2*y(1) - 16*s^2*y(3) - 32*s^2*y(4) + 32*s^2*y(5) + 32*r*s*y(1) + 32*r*s*y(2) - 64*r*s*y(4);
                       -(4*r + 4*s - 3)*(x(2) - x(3) - 4*x(4) + 4*x(6) - 4*r*x(2) + 4*r*x(4) + 4*r*x(5) - 4*r*x(6) + 4*s*x(3) + 4*s*x(4) - 4*s*x(5) - 4*s*x(6)),  (4*r - 1)*(3*x(1) + x(3) - 4*x(6) - 4*r*x(1) + 4*r*x(4) - 4*r*x(5) + 4*r*x(6) - 4*s*x(1) - 4*s*x(3) + 8*s*x(6)), -(4*s - 1)*(3*x(1) + x(2) - 4*x(4) - 4*r*x(1) - 4*r*x(2) + 8*r*x(4) - 4*s*x(1) + 4*s*x(4) - 4*s*x(5) + 4*s*x(6)), 12*x(1) + 4*x(3) - 16*x(6) - 28*r*x(1) + 4*r*x(2) - 8*r*x(3) - 16*r*x(5) + 48*r*x(6) - 28*s*x(1) - 20*s*x(3) + 48*s*x(6) + 16*r^2*x(1) - 16*r^2*x(2) + 32*r^2*x(5) - 32*r^2*x(6) + 16*s^2*x(1) + 16*s^2*x(3) - 32*s^2*x(6) + 32*r*s*x(1) + 32*r*s*x(3) - 64*r*s*x(6), 16*r*x(4) - 4*r*x(2) - 12*r*x(1) + 12*s*x(1) + 4*s*x(3) - 16*s*x(6) + 16*r^2*x(1) + 16*r^2*x(2) - 32*r^2*x(4) - 16*s^2*x(1) - 16*s^2*x(3) + 32*s^2*x(6), 16*x(4) - 4*x(2) - 12*x(1) + 28*r*x(1) + 20*r*x(2) - 48*r*x(4) + 28*s*x(1) + 8*s*x(2) - 4*s*x(3) - 48*s*x(4) + 16*s*x(5) - 16*r^2*x(1) - 16*r^2*x(2) + 32*r^2*x(4) - 16*s^2*x(1) + 16*s^2*x(3) + 32*s^2*x(4) - 32*s^2*x(5) - 32*r*s*x(1) - 32*r*s*x(2) + 64*r*s*x(4)];
        
        B = [B_SF(1,1),0,B_SF(1,2),0,B_SF(1,3),0,B_SF(1,4),0,B_SF(1,5),0,B_SF(1,6),0;
             0, B_SF(2,1),0,B_SF(2,2),0,B_SF(2,3),0,B_SF(2,4),0,B_SF(2,5),0,B_SF(2,6);
             B_SF(1,1),B_SF(2,1),B_SF(1,2),B_SF(2,2),B_SF(1,3),B_SF(2,3),B_SF(1,4),B_SF(2,4),B_SF(1,5),B_SF(2,5),B_SF(1,6),B_SF(2,6)];
end