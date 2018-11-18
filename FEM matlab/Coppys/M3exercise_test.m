% FE program for MEFEM class: 2D plane stress Q4 elements
% MMP, rev. 4, 26/9-2018.


%%%%%%%%%%%%%%%%%%%% define model %%%%%%%%%%%%%%%%%%%%%%%%%
W  = 30;        % model width
H  = 30;        % model height
h  = 1;         % model thickness

nw = 4;         % no. elements in x-direction
nh = 4;         % no. elements in y-direction
nGP = 2;        % no. Gauss points

E  = 210000;    % Youngs modulus
nu = 0.25;       % Poissons ratio

% generate mesh (node and element tables)
[nodes,elems,ndof] = Q4mesh(W,H,nw,nh);

% set boundary conditions
BC = false(ndof,1);
% fixed dofs (at x = 0 and y = 0), all others are free
BC(find(nodes(:,1) == 0)*2-1) = 1; 
BC(find(nodes(:,2) == 0)*2) = 1;
% apply nodal loads (forces directly on nodes)
P = zeros(ndof,1);
P([10,20,30,40,50]) = -1000;
R = P;


%%%%%%%%%%%%%%%%%%%%%% plot model %%%%%%%%%%%%%%%%%%%%%%%%%%
% (comment out what you don't need)
init_plot(nodes); % always keep this one!
plot_elements(elems,nodes)
plot_elem_numbers(elems,nodes)
% plot_elem_cs(elems,nodes)
plot_nodes(nodes)
plot_node_numbers(nodes)
plot_dof_numbers(nodes)
plot_BCs(nodes,BC)
plot_loads(nodes,P)



%%%%%%%%%%%%%%%%%%%%%% setup & solve %%%%%%%%%%%%%%%%%%%%%%%
% constitutive matrix for plane stress, Cook eq. (3.1-3)
Em = [ 1 nu 0;
      nu  1 0;
       0  0 (1-nu)/2]*E/(1-nu^2);

% assembly of global stiffness matrix
K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP)
 
% solve for displacements and reaction
[D,R] = solve(K,R,BC);



%%%%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%
print_nodal_results(D,R,BC) % show results command window
plot_reaction_forces(nodes,D,R)
plot_displacement(elems,nodes,D*100,'mesh') %last argument: 'UX', 'UY', 'U' or 'mesh'





%%%%%%%%%%%%%%%%%% FE functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = kQ4(a,b,Em,h,nGP)
% Determines the element stiffness matrix for a bilinear 
% rectangular Q4 element using numerical integration. 
% INPUTS:
%   a   =    half element width
%   b   =    half element height
%   Em  =    constitutive matrix [3x3]
%   h   =    element thickness
%   nGP =    number of Gauss points for numerical integration
% OUTPUTS:
%   k   =    the 8x8 stiffness matrix

k = zeros(8);
[GP,W] = Gauss(nGP);

for i = 1:nGP
    J = a*b;
    for j = 1:nGP
        
        B = Bmat(a*GP(i),b*GP(j),a,b);       

        k = k + B'*Em*B*h*J*W(i)*W(j);        
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

function B = Bmat(x,y,a,b)
% Evaluates the 3x8 strain-displacement matrix [B] for a bilinear 
% rectangular element Q4. Cook eq. (3.6-6)
%
% INPUTS:
%   x =    x-coordinate at which to evaluate B
%   y =    y-coordinate at which to evaluate B
%   a =    half element width
%   b =    half element height
%
% OUTPUTS: 
%   B =    strain displacement matrix 

B = 1/(4*a*b) * [-(b-y) 0   (b-y)   0   (b+y)   0   -(b+y)   0;
                0   -(a-x)  0   -(a+x)  0   (a+x)   0   (a-x);
                -(a-x)  -(b-y)  -(a+x)  (b-y)   (a+x)   (b+y)   (a-x)   -(b+y)];
        
end

%%
function [D,R] = solve(K,R,BC)
% solve FE problem, Cook eq. (2.7-3, 2.7-2b)
% D = vector of displacements (both free and fixed)
% R = vector of nodal forces (applied and reactions)

% Splitting global stiffness matrix
K11 = K(~BC,~BC);
K12 = K(~BC,BC);
K21 = K(BC,~BC);
K22 = K(BC,BC);

% Defining the constrained deformation and forces
Dc = zeros(sum(BC),1);
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

%%
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
        
         x=nodes(nod,1);
         y=nodes(nod,2);
         
         a = x(2)-x(1);
         b = y(3)-y(2);
         
        dof_x = nod*2-1;
        dof_y = nod*2;
        
        dof = reshape([dof_x;dof_y],1,[]);
        
        k = kQ4(a,b,Em,h,nGP);
        

        K(dof,dof) = K(dof,dof) + k;

    end
        
end





% Meshing: all taken care of. you're welcome.
function [nodes,elems,ndof] = Q4mesh(W,H,nw,nh)

% Square panel mesh generation.
%
% INPUTS:
%  W  = width of model
%  H  = height of model
%  nw = no. of elements in x-direction
%  nh = no. of elements in y-direction
%
% OUTPUTS:
% nodes = node table: 1=x-coordinate, 2=y-coordinate, (row no. = node no.)
% elems = element table: 1=node i, 2=node j, ... (row no. = elem no.) 

    % no. of dof's
    ndof = 2*(nw+1)*(nh+1);

    % no. of nodes
    nn = ndof/2;

    % no. elements
    ne = nw*nh;

    % initialize arrays
    nodes = zeros(nn,2);
    elems = zeros(ne,4);

    % grid of nodal coordinates
    x_grid = linspace(0,W,nw+1);
    y_grid = linspace(0,H,nh+1);

    % setup node table
    n = 0; % node number
    for i = 1:(nw+1)
        for j=1:(nh+1)
            n = n+1;
            nodes(n,1) = x_grid(i);
            nodes(n,2) = y_grid(j);
        end
    end

    % setup element table
    e_no = 0;
    for j=1:nh
        for i = 1:nw
            e_no = e_no+1;

            n1 = (nh+1)*i + j -(nh+1);
            n2 = (nh+1)*i + j;
            n3 = n2 + 1;
            n4 = n1 + 1;

            elems(e_no,1:4) = [n1 n2 n3 n4];
        end
    end

end








%% %%%%%%%%%%%%%%%% plotting functions %%%%%%%%%%%%%%%%%%%%%
% don't worry about these. 
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

function plot_reaction_forces(nodes,D,R)

    for dof = 1:length(D)
        if D(dof)==0 % fixed dof
            
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
