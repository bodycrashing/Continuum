% FE program for MEFEM class: 2D plane stress Q4 elements
% MMP, rev. 4, 26/9-2018.
clc; clear all; close all;

%% %%%%%%%%%%%%%%%%%% define model %%%%%%%%%%%%%%%%%%%%%%%%%
W  = 30;        % Unit zell widht and height [mm]
R  = 20;        % Radius of hole [mm]
h  = 1;         % model thickness

nc = 4;         % no. elements in on hole curve
no = 4;         % no. elements in allong diagonal of hole curve
% nc MUST!! be even numbers
nGP = 2;        % no. Gauss points

E  = 2.1e11;    % Youngs modulus
nu = 0.25;       % Poissons ratio

%% %%%%%%%% Generate Mesh, Constraints and Loads %%%%%%%%%%%%%%%%%%%%%%%

% generate mesh (node and element tables)
[nodes,elems,ndof] = Q4mesh_UnitCell(W,R,nc,no);

% set boundary conditions
BC = false(ndof,1);
MC = false(ndof,1);
% fixed dofs (at x = 0 and y = 0), all others are free
BC(find(nodes(:,1) == 0)*2-1) = 1;
BC(find(nodes(:,2) == 0)*2) = 1;
%BC(find(nodes(:,1) == W)*2-1) = 1;
BC(find(nodes(:,2) == W)*2) = 1;
% Forced siplacements
Dc = zeros(ndof,1);
Dc(find(nodes(:,2) == W)*2) = 2;

% Line load at top of unit cell. Value in N/mm2
poly = [50];

% Calculate nodal lodas from line load
%P = element_loads(poly,elems,nodes,h,W,nGP);
P = zeros(ndof,1);
R = P;
% Multivariable constraints
MC(find(nodes(:,1) == W)*2-1) = 1;
C = mulit_constraints(elems,nodes,MC);


%% %%%%%%%%%%%%%%%%%%%% setup & solve %%%%%%%%%%%%%%%%%%%%%%%
% constitutive matrix for plane stress, Cook eq. (3.1-3)
Em = [ 1 nu 0;
    nu  1 0;
    0  0 (1-nu)/2]*E/(1-nu^2);

% Assembly of global stiffness matrix
K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP);

% solve for displacements and reaction
[D,R] = solve(K,R,BC,Dc,C);

% Calculating stresses in the elements
[eps,sigma_val] = StressCal(elems,nodes,D,E);


%% %%%%%%%%%%%%%%%%%%%% plot model %%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Plot
% (comment out what you don't need)
init_plot(nodes); % always keep this one!
plot_elements(elems,nodes)
plot_elem_numbers(elems,nodes)
%plot_elem_cs(elems,nodes)
plot_nodes(nodes)
plot_node_numbers(nodes)
plot_dof_numbers(nodes)
plot_BCs(nodes,BC)
%plot_loads(nodes,P)
%plot_dist_load(elems,nodes,poly,W)

% Results Plot
%print_nodal_results(D,R,BC) % show results command window
plot_reaction_forces(nodes,D,R,BC)


plot_stress(elems,nodes,sigma_val)
figure()
plot_displacement(elems,nodes,D,'mesh') %last argument: 'UX', 'UY', 'U' or 'mesh'


%% %%%%%%%%%%%%%%%% FE functions %%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = kQ4Iso(x_vec,y_vec,Em,h,nGP)
% Determines the element stiffness matrix for a isoparametric
% Q4 element using numerical integration.
% INPUTS:
%   x_vec   =    vector containig x coordinates of element nodes
%   y_vec   =    vector containig y coordinates of element nodes
%   Em  =    constitutive matrix [3x3]
%   h   =    element thickness
%   nGP =    number of Gauss points for numerical integration
% OUTPUTS:
%   k   =    the 8x8 stiffness matrix

% Initialize k
k = zeros(8);
% Get Gauss points and weights
[GP,W] = Gauss(nGP);

for i = 1:nGP    
    for j = 1:nGP
        
        % xi and eta are set to be the gaus points
        xi = GP(i);
        eta = GP(j);
        % Evaluating the B matrix at xi and eta       
        [B,J_det] = Bmat(xi,eta,x_vec,y_vec);
        
        k = k + B'*Em*B*h*J_det*W(i)*W(j);
    end
end
end

function [GP,W] = Gauss(nGP)
% Returns sampling points and weights for Gauss quadrature of order nGP
% Cook, Table 6.3-1
% INPUTS:
%   nGP =    order of quadrature
% OUTPUTS:
%   GP  =    sampling points
%   W   =    weight functions


w_l = [2,0,0;
    1, 1,0;
    5/9, 8/9, 5/9];

GP_l = [0,0,0;
    -1/sqrt(3),1/sqrt(3),0;
    -sqrt(3/5),0,sqrt(3/5)];

GP = zeros(nGP,1);

GP(:) = GP_l(nGP,1:nGP);
W = w_l(nGP,1:nGP);



end

function N = Q4shape(xi,eta)
% evaluates the shape functions at coordinates x,y
% a,b    - half-length of element
% x,y    - coordinates of sampling point
% N      - shape functions, size 2x8

% This function is the same as for regular Q4, but for isoparametric Q4
% a and b are always equal to 1. x and y are replaced by xi and eta;
a = 1;
b = 1;

N1 = 1/(4*a*b)*(a-xi)*(b-eta);
N2 = 1/(4*a*b)*(a+xi)*(b-eta);
N3 = 1/(4*a*b)*(a+xi)*(b+eta);
N4 = 1/(4*a*b)*(a-xi)*(b+eta);

N = kron([N1 N2 N3 N4],eye(2));

end

function [B,J_det] = Bmat(xi,eta,x_vec,y_vec)
% Calculating B matrix for isoparametric Q4 following principles from 
% Cook p. 205-208.
% xi and eta = points to evaluate B at. Must be between -1 and 1
% x_vec and y_vec = Vertical vectors containing coordinates of element nodes

% Calculating the Jacobian matrix, eq 6.2-6 in cook
J1 = 1/4 * [-(1-eta) (1-eta) (1+eta) -(1+eta);
    -(1-xi) -(1+xi)  (1+xi)   (1-xi)];

J = J1*cat(2,x_vec,y_vec);
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

end

function [D,R] = solve(K,R,BC,Dc,C)
% solve FE problem, Cook eq. (2.7-3, 2.7-2b)
% D = vector of displacements (both free and fixed)
% R = vector of nodal forces (applied and reactions)

O = zeros(size(C),1);

% Splitting global stiffness matrix
K11 = K(~BC,~BC);
K12 = K(~BC,BC);
K21 = K(BC,~BC);
K22 = K(BC,BC);

% Defining the constrained deformation and forces
Dc = Dc(BC);
Rc = R(~BC);

% Calculating free deformations and forces
Dx = K11\(Rc-K12*Dc);

Rx = K21*Dx + K22*Dc;

% Filling D and R vector
D = zeros(length(BC),1);
R = zeros(length(BC),1);

D(~BC)=Dx;
D(BC) = Dc;

R(~BC) = Rc;
R(BC) = Rx;


end

function K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP)
% Assembles the global/structure stiffeness matrix K from element stiffness matrices
%
% INPUTS:
%   nodes = node table: 1=x-coordinate, 2=y-coordinate, (row no. = node no.)
%   elems = element table: 1=node i, 2=node j, ... (row no. = elem no.)
%   Em    = constitutive matrix
%   h     = thickness
%   ndof  = total number of dofs in model
%   nGP   = no. Gauss points for integration
%
% OUTPUTS:
%   K     = global stiffness matrix [ndof x ndof]


% develop this function from scratch

K = zeros(ndof,ndof);
% K = spalloc(ndof,ndof,max_nonzero_entries) % this is better for large models.

for i = 1:size(elems,1)
    nod = elems(i,:);
    
    x_vec=nodes(nod,1);
    y_vec=nodes(nod,2);
    
    
    dof_x = nod*2-1;
    dof_y = nod*2;
    
    dof = reshape([dof_x;dof_y],1,[]);
    
    
    k = kQ4Iso(x_vec,y_vec,Em,h,nGP);
    
    
    K(dof,dof) = K(dof,dof) + k;
    
end

end

function [eps,sigma_val] = StressCal(elems,nodes,D,E)
% Calculating stresses in the elements
% D = Displacements
% E = Modulus of elasticity
% Output:
% eps = 3rd order tensor containing strains. 
%   First index = strain type [e_x, e_y, e_xy] . Values from 1 to 3
%   Second index = node number in element. 1 to 4
%   Third index = element number. 1 to number of elements

for i = 1:size(elems,1) % Looping over number of elements
    
    % Calculating dof belonging to current element
    nod = elems(i,:);
    dof_x = nod*2-1;
    dof_y = nod*2;    
    dof = reshape([dof_x;dof_y],1,[]);
    % Fetching displacements belonging to current element
    d = D(dof);
    % Getting x and y coordinates for element nodes
    x_vec=nodes(nod,1);
    y_vec=nodes(nod,2); 
    
    xi = [-1 1 1 -1];
    eta = [-1 -1 1 1];
    % Looping over number of nodes pr element
    for j = 1:size(elems,2) 
        % Evaluating B at node position  
        B =  Bmat(xi(j),eta(j),x_vec,y_vec);
        % Calculating stresses and strains
        eps(:,j,i) = B*d;
        sigma_val(:,j,i) = E*eps(:,j,i);
        
    end
end

end

function R = element_loads(poly,elems,nodes,h,W,nGP)
% Calculate the nodal loads from the distributed loads on the element
% INPUTS:
% poly  = polynomial coefficients of applied stress function in y-dir,
% elems = element table (one Q4 element only)
% nodes = 
% h     = thickness of element
% W     = Width of unit cell
% nGP   = no. Gauss points
% OUTPUTS:
% R    = vector of nodal loads

% init load vector
R = zeros(size(nodes,1)*2,1);

% Looping over elements
for i = 1:size(elems(:,1)) 
    
    % Only If element is located at upper boundary
    % nodal forces should be calculated for the element
    if any(sum(nodes(elems(i,:),2)==W)==2)
       
        % Initialize nodal load for element
        re = zeros(size(elems,2)*2,1);
        % Getting x and y coordinates of element nodes
        x_vec = nodes(elems(i,:),1);
        y_vec = nodes(elems(i,:),2);       
        c = reshape([x_vec';y_vec'],1,[])';
        % Creating function expression for distributed load       
        phi = @(x) [0, polyval(poly,x)/h]';
        
        % The scaling factor J is equal to the length of top 
        % side width of the element
        top_nodes = nodes(elems(i,:),2)==W;
        top_nodes_x_val = nodes(elems(i,top_nodes),1);      
        J = abs(diff(top_nodes_x_val))/2;
        
        [GP,Wg] = Gauss(nGP);
        
        % Using Gauss Quadrature to calculate loads
        for j = 1:nGP
            % Inside the nodes allong the top of the unit cell 
            % xi is equal to 1 and eta varies with x.
            eta = GP(j);
            xi = 1;
            % Evaluating the shape function at the Gauss point
            N = Q4shape(xi,eta);
            % Getting the x and y coordinates of Gauss point
            xy = N*c;
            % The nodal load is calculated
            re = re + N'*phi(xy(1)) *h*J* Wg(j);
            
        end  
        % The element load are inserted in to the global load vector R
        nod = elems(i,:);
        dof_x = nod*2-1;
        dof_y = nod*2;
        dof = reshape([dof_x;dof_y],1,[]);
        
        R(dof) = R(dof) + re;
    end
end

end

function C = mulit_constraints(elems,nodes,MC)

test = [1 0 -1];
ndof = length(nodes)*2;
C = zeros(sum(MC)-1,ndof);

for i = 1:sum(MC)-1
    
    C(i,(i*2-1:i*2+1)) = test;
    
end


end

%% %%%%%%%%%%%%%%%%%% plotting functions %%%%%%%%%%%%%%%%%%%%%

function init_plot(nodes)
% setup figure window and plot undeformed model

figure('name','MEFEMquad','NumberTitle','off');
centerfig;

% initialize plot
plot(0,0,'+k'); % origin
hold on
axis equal
xlabel('x');
ylabel('y');


% adjust plot limits
nx = nodes(:,1);
ny = nodes(:,2);
rx = max(nx)-min(nx);
ry = max(ny)-min(ny);

c = max(5,max(max(nx),max(ny))/5);

if rx>ry
    xlim([min(nx)-c max(nx)+c])
else
    ylim([min(ny)-c max(ny)+c])
end

end

function print_nodal_results(D,R,BC)

% print nodal results to command window
fprintf('NODAL RESULTS:\n')
fprintf('node  x-dof \t  D_x \t\t   R_x\t\ty-dof \t  D_y \t\t   R_y\n')
fprintf('----------------------------------------------------------------------\n')

nn = size(R,1)/2;

for n = 1:nn
    
    % dof no's
    x = 2*n-1;
    y = 2*n;
    
    % put '=' in front of fixed displacements
    if BC(x)
        Dcx = '=';
    else
        Dcx = ' ';
    end
    if BC(y)
        Dcy = '=';
    else
        Dcy = ' ';
    end
    
    fprintf('%4i  %5i \t%s% 6.2e \t % 6.2e \t %4i \t%s% 6.2e \t % 6.2e \n',n,x,Dcx,D(x),R(x),y,Dcy,D(y),R(y))
end

end

function plot_BCs(nodes,BC)
% plot boundary conditions

for dof = 1:length(BC)
    if BC(dof)==1 % fixed dof
        
        % associated node
        n = ceil(dof/2);
        nvec = nodes(n,1:2)';
        
        % draw triangle
        if rem(dof,2)==1 % x-dir
            v_verts = [nvec';
                nvec' + 3*[-0.5  0.25];
                nvec' + 3*[-0.5 -0.25]];
            patch('faces',[1 2 3],'vertices',v_verts,'edgecolor','k','facecolor','none','tag','Boundary conditions');
        else % y-dof
            h_verts = [nvec';
                nvec' + 3*[0.25 -0.5];
                nvec' + 3*[-0.25 -0.5]];
            patch('faces',[1 2 3],'vertices',h_verts,'edgecolor','k','facecolor','none','tag','Boundary conditions');
        end
    end
end

end

function plot_loads(nodes,P)
% nodal loads, P

for dof = 1:size(P,1)
    if P(dof)~=0
        
        % associated node
        n  = ceil(dof/2);
        nvec = nodes(n,1:2)';
        
        % normalized load
        Pn = 2*P(dof)/max(abs(P));
        
        % vector
        if rem(dof,2)==1 % x-dir
            drawArrow([nvec(1) nvec(1)+Pn],[nvec(2) nvec(2)],'r','linewidth',2,'tag','Loads');
        else
            drawArrow([nvec(1) nvec(1)],[nvec(2) nvec(2)+Pn],'r','linewidth',2,'tag','Loads');
        end
    end
end

end

function plot_reaction_forces(nodes,D,R,BC)

for dof = 1:length(D)
    if BC(dof) % fixed dof
        
        % associated node
        n = ceil(dof/2);
        node_vec = nodes(n,1:2)';
        
        % normalized reaction
        Rn = 2*R(dof)/max(abs(R));
        
        % normalized reaction vector
        if rem(dof,2)==1 % x-dir
            drawArrow([node_vec(1) node_vec(1)+Rn],[node_vec(2) node_vec(2)],'m','linewidth',2,'tag','Reactions');
        else
            drawArrow([node_vec(1) node_vec(1)],[node_vec(2) node_vec(2)+Rn],'m','linewidth',2,'tag','Reactions');
        end
        
    end
end

end

function plot_nodes(nodes)
% nodes only

for i = 1:size(nodes,1)
    plot(nodes(i,1),nodes(i,2),'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor',[0.5 0.5 0.5]);
end

end

function plot_node_numbers(nodes)
% node numbers

for i = 1:size(nodes,1)
    text(nodes(i,1)+1,nodes(i,2)+1,num2str(i),'color','r');
end

end

function plot_dof_numbers(nodes)

for i = 1:size(nodes,1)
    text(nodes(i,1)+1,nodes(i,2)-1,['(' num2str(i*2-1) ',' num2str(i*2) ')'],'color','b','fontsize',7);
end

end

function plot_elements(elems,nodes)
% plot elements

for i = 1:size(elems,1)
    
    % nodes on element
    n = elems(i,1:4);
    
    % element retangles
    verts = nodes(n,:);
    patch('faces',1:4,'vertices',verts,'facecolor','c');
    
end

end

function plot_elem_cs(elems,nodes)
% plot element numbers + CS

for i = 1:size(elems,1)
    
    % nodes on element
    n = elems(i,1:4);
    
    % element CS & number
    mid = [0 0]';
    for j = 1:4
        vec = [nodes(n(j),1) nodes(n(j),2)]';  % vectors to nodes
        mid = mid + vec;
    end
    mid = mid/4;
    
    drawArrow([mid(1) mid(1)+1],[mid(2) mid(2)],'b');
    drawArrow([mid(1) mid(1)]  ,[mid(2) mid(2)+1],'b');
    
end

end

function plot_elem_numbers(elems,nodes)
% plot element numbers + CS

for i = 1:size(elems,1)
    
    % nodes on element
    n = elems(i,1:4);
    
    % element CS & number
    mid = [0 0]';
    for j = 1:4
        vec = [nodes(n(j),1) nodes(n(j),2)]';  % vectors to nodes
        mid = mid + vec;
    end
    mid = mid/4;
    
    text(mid(1),mid(2),num2str(i),'horizontalalignment','center');
    
end

end

function plot_displacement(elems,nodes,D,plot_mode)

% x/y dof numbers
nn = length(nodes);
x_dof = 1:2:nn*2;
y_dof = 2:2:nn*2;

% pick out displacements for x- and y-dir, respectively
UX = D(x_dof);
UY = D(y_dof);
U  = sqrt(UX.^2 + UY.^2);

% add displacements to node coordinates
nodes = nodes + [UX UY];

switch plot_mode
    case 'UX'
        phi = UX;
    case 'UY'
        phi = UY;
    case 'U'
        phi = U;
end

if strcmp(plot_mode,'mesh')
    patch('faces',elems,'vertices',nodes,'facecolor','none');
else
    patch('faces',elems,'vertices',nodes,'facecolor','interp','FaceVertexCData',phi);
    cb = colorbar;
    colormap(jet(8));
    cb.Title.String = plot_mode;
end

end

function v_rot = rot(v,theta)
% Rotates v by the angle theta

A = [cos(theta) -sin(theta);
    sin(theta)  cos(theta)];

v_rot = A*v;

end

function h = drawArrow(x,y,varargin)
% define homemade arrow function

% unit vector in arrow direction
u(1,1) = x(2)-x(1);
u(2,1) = y(2)-y(1);
u = u/norm(u);

% arrow body
h(1) = plot([x(1) x(2)],[y(1) y(2)],varargin{:});

% head at arrow end (2)
a = [x(2) y(2)]';
b = a-0.1*u;
c = a-b;
d = b + rot(c,pi/2);
e = b - rot(c,pi/2);
h(2) = plot([d(1) a(1) e(1)],[d(2) a(2) e(2)],varargin{:});

end

function plot_dist_load(elems,nodes,poly,W)
    % distributed load on top surface of element 1
    
    for k = 1:size(elems,1)
    
    if any(sum(nodes(elems(k,:),2)==W)==2)
    ei = k;
        
    n4 = elems(ei,3);
    n3 = elems(ei,2);
    
    n4_vec = nodes(n4,1:2)';
    n3_vec = nodes(n3,1:2)';
    
    % half side length
    a = abs( n4_vec(1) - n3_vec(1) )/2;
    
    % vector along top surface of element
    e = (n3_vec-n4_vec); 
    u = e/norm(e); % unit

    % distributed load
    na = 10;
    x = linspace(n4_vec(1),n3_vec(1),na+1);
    q = polyval(poly,x);
    
    
    % draw arrows over element length
    for j = 0:na
        p_beg = n4_vec + (j/na)*e; % arrow end point
        p_end = p_beg - 2*rot(u,pi/2)*q(j+1)/max(abs(q)); % arrow begin point
        drawArrow([p_beg(1) p_end(1)],[p_beg(2) p_end(2)],'b','linewidth',1.5);
    end

    % add line above arrows
    y = n4_vec(2);
    plot(x,y+2*q/max(q),'b','linewidth',1);
    end
    end
end

function plot_stress(elems,nodes,sigma_val)
% Plotting stresses 
%figure()
figure()
ax = axes();
axis equal
hold on
colorbar

for i = 1:size(elems,1)
    for j = 1:size(elems,2)
        % Calculating von mises stress
        sigma_vm(j,i) = sqrt( sigma_val(1,j,i)^2 + sigma_val(2,j,i)^2 - sigma_val(1,j,i)*sigma_val(2,j,i) + 3*sigma_val(3,j,i));
        
    end
    eface = 1:size(elems,2);
    n = elems(i,eface);
    enodes = nodes(n,:);
    hold on
    
    if ~isreal(sigma_vm)
       sigma_vm = real(sigma_vm)
    end
    
    h(i) = patch('faces',eface,'vertices',enodes,...
          'facecolor','interp','FaceVertexCData',sigma_vm(:,i));


    
end
t1 = hgtransform('Parent',ax);
t2 = hgtransform('Parent',ax);
t3 = hgtransform('Parent',ax);
t4 = hgtransform('Parent',ax);

set(h,'Parent',t1)
h2 = copyobj(h,t2);
h3 = copyobj(h,t3);
h4 = copyobj(h,t4);

Rx = makehgtform('xrotate',pi);
Ry = makehgtform('yrotate',pi);

set(t2,'Matrix',Rx)
set(t3,'Matrix',Rx*Ry)
set(t4,'Matrix',Ry)





end


