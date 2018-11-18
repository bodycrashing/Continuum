% FE program for MEFEM class: 2D plane stress Q4 elements
% MMP, rev. 3, 31/7-2017.
clc; close all; clear all

%%%%%%%%%%%%%%%%%%%% define model %%%%%%%%%%%%%%%%%%%%%%%%%
W  = 50;        % model width
H  = 10;        % model height
h  = 1;         % model thickness

nw = 5;         % no. elements in x-direction
nh = 2;         % no. elements in y-direction
nGP = 2;        % no. Gauss points

E  = 210000;    % Youngs modulus
nu = 0.3;       % Poissons ratio

% generate mesh (node and element tables)
[nodes,elems,ndof] = Q4mesh(W,H,nw,nh);

% set boundary conditions
BC = false(ndof,1);
BC([1 2 3 5],1) = 1; % fixed dofs, all others are free

% apply nodal loads (forces directly on nodes)
P = zeros(ndof,1);
P(end) = -5000;

R = P;


%%%%%%%%%%%%%%%%%%%%%% plot model %%%%%%%%%%%%%%%%%%%%%%%%%%
% (comment out, what you don't need)
init_plot(nodes); % always keep this one!
plot_elements(elems,nodes)
%plot_elem_numbers(elems,nodes)
%plot_elem_cs(elems,nodes)
plot_nodes(nodes)
%plot_node_numbers(nodes)
plot_dof_numbers(nodes)
plot_BCs(nodes,BC)
plot_loads(nodes,P)



%%%%%%%%%%%%%%%%%%%%%% setup & solve %%%%%%%%%%%%%%%%%%%%%%%
% constitutive matrix for plane stress
Em = [ 1 nu 0;
      nu  1 0;
       0  0 (1-nu)/2]*E/(1-nu^2);

% assembly of global stiffness matrix
K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP);

% solve for displacements and reaction
[D,R] = solve(K,R,BC);

% calculate stress and strains etc.
res = post_process(elems,nodes,D,Em);

%%%%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%
print_nodal_results(D,R,BC) % show results command window
plot_reaction_forces(nodes,D,R)
plot_displacement(elems,nodes,D,'mesh') %last argument: 'UX', 'UY', 'U' or 'mesh'
plot_results(elems,nodes,res,'SX')



%%%%%%%%%%%%%%%%%% FE functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [GP,W] = Gauss(n)
% Returns sampling points and weights for for Gauss quadrature of order n
% Cook, Table 6.3-1
%    xi:      sampling points
%     W:      weight functions
%     n:      order of quadrature

    switch n 
        case 1
            GP = 0;
            W  = 2;
        case 2
            GP = [-sqrt(1/3) sqrt(1/3)];
            W  = [1 1];
        case 3
            GP = [-sqrt(0.6) 0 sqrt(0.6)];
            W  = [5/9 8/9 5/9];
    end

end

function k = kQ4(x,y,Em,h,nGP)
% Function for calculating the local stiffness matrix
% for a bilinear rectangular Q4 element using numerical
% integration. 
%   k:      the 8x8 stiffness matrix
%   x:      global element coordinate x (dimension 4x1).
%   y:      global element coordinate y (dimension 4x1).
%   Em:     3x3 constitutive matrix
%   h:      element thickness
%   nGP:    number of Gauss points for numerical integration

    % get sampling points and weights for numerical integration
    [GP,W] = Gauss(nGP);

    % sampling points
    xi = GP;
    eta = GP;

    % pre-allocating empty 8x8 matrix
    k = zeros(8,8);

    for i = 1:nGP
        for j = 1:nGP
            % Evaluating B matrix at sampling point (i,j)
            [B,J] = Bmat(x,y,xi(i),eta(j));
            
            % Calculating contribution from current sampling point
            dk = B'*Em*B*h*J*W(i)*W(j);

            % Summing up contribution with the stiffness matrix
            k = k + dk;        
        end
    end

end

function [B,J] = Bmat(x,y,xi,eta)
% Evaluates the 3x8 strain-displacment matrix [B] for a bilinear 
% rectangular element Q4. Cook eq. (3.6-6)
%   B:      strain displacement matrix
%   J:      Jacobian
%   x:      global element coordinate x (dimension 4x1).
%   y:      global element coordinate y (dimension 4x1).
%   xi:     local isoparametric coordinate xi (scalar).
%   eta:    local isoparametric coordinate eta (scalar).
    
    B1 = [1 0 0 0;
          0 0 0 1;
          0 1 1 0];                         % Plane strain-displacement relations eq. 6.2-9
    
    % Shape functions are given by eq 6.2-3
    N1xi = -1/4*(1-eta);                   % Partial derivative of N1 with respect xi
    N2xi =  1/4*(1-eta);                   % Partial derivative of N2 with respect xi
    N3xi =  1/4*(1+eta);                   % Partial derivative of N3 with respect xi
    N4xi = -1/4*(1+eta);                   % Partial derivative of N4 with respect xi
    
    N1eta = -1/4*(1-xi);                   % Partial derivative of N1 with respect eta
    N2eta = -1/4*(1+xi);                   % Partial derivative of N2 with respect eta
    N3eta =  1/4*(1+xi);                   % Partial derivative of N3 with respect eta
    N4eta =  1/4*(1-xi);                   % Partial derivative of N4 with respect eta
      
    Jmatrix = [N1xi N2xi N3xi N4xi;
               N1eta N2eta N3eta N4eta]*...
                  [x(1) y(1);
                   x(2) y(2);
                   x(3) y(3);
                   x(4) y(4)];              % Jacobian matrix eq. 6.2-5
               
    J = det(Jmatrix);                       % Jacobian eq. 6.2-8
    
    Gamma = inv(Jmatrix);                   % Inverse of Jacobian matrix eq. 6.2-7
      
    B2 = zeros(4);
    B2(1:2,1:2) = Gamma;
    B2(3:4,3:4) = Gamma;                    % Expanded version of eq. 6.2-7 due to the two displacement fields u,v

    B3 = [N1xi 0 N2xi 0 N3xi 0 N4xi 0;
          N1eta 0 N2eta 0 N3eta 0 N4eta 0;
          0 N1xi 0 N2xi 0 N3xi 0 N4xi;
          0 N1eta 0 N2eta 0 N3eta 0 N4eta]; % Local xi-eta strain displacement relations

    B = B1*B2*B3;                           % Final strain displacement matrix eq. 6.2-9 to 6.2-11
end

function K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP)
% Assemble global stiffness matrix by summation over elements
%
%   INPUTS:
%   elems = element table
%   Em    = constitutive matrix
%   h     = element thickness
%   ndof  = no. dof's
%   nGP   = no. Gauss points
%
%   OUTPUTS:
%   K     = global structure stiffness matrix

    % initialize
    ne = size(elems,1);   % no. elements
    K  = zeros(ndof);     % empty global stiffness matrix
    
    for i = 1:ne

        % node numbers of current element
        n = elems(i,1:4);

        % dof numbers
        dof = [2*n(1)-1 2*n(1) 2*n(2)-1 2*n(2) 2*n(3)-1 2*n(3) 2*n(4)-1 2*n(4)]';
        
        % n coordinates
        x = nodes(n,1);
        y = nodes(n,2);
        
        % element stiffness matrix in local cs, Cook (2.3-5/6)
        k = kQ4(x,y,Em,h,nGP);

        % assemble structure stiffness matrix, Cook (2.5-4)
        K(dof,dof) = K(dof,dof) + k;
        
    end
        
end

function [D,R] = solve(K,R,BC)
% solve FE problem
% D = vector of displacements (both free and fixed)
% R = vector of nodal forces (applied and reactions)

    % initialize outputs
    ndof = size(K,1);
    D    = zeros(ndof,1);
    
    % indices of fixed and free dofs:
    fixed = BC;
    free  = ~BC;
    
    % split up K & R
    K11 = K(free,free);
    K21 = K(fixed,free);
    K12 = K(free,fixed);
    K22 = K(fixed,fixed);
    Rc  = R(free);
    
    % solve for displacements (Dc=0)
    Dc = D(BC);
    Dx = K11\(Rc-K12*Dc);

    % solve for reactions
    Rx  = K21*Dx + K22*Dc;

    % put results into array at correct dof places
    D(free)  = Dx;
    R(free)  = Rc;
    R(fixed) = Rx;
   
end

function res = post_process(elems,nodes,D,Em)

    % numbers of elements and nodes
    ne = size(elems,1);
    nn = size(nodes,1);

    % loop over elements: calculate nodal results for each element (unaveraged)
    for i = 1:ne
        
        
        
        
        
        
        
        % put code here...
        
        
        
        
        
        
        
        
        
        
        % loop over element nodes
        for j = 1:4
            

            
            
            
            
            % replace these 2 lines with your calculations
            eps   = zeros(3,1);
            sigma = zeros(3,1);
                        
            
            
            
            
            
            
            
            
            
            % store results per element
            res.element(i).EX(j)  = eps(1);
            res.element(i).EY(j)  = eps(2);
            res.element(i).EXY(j) = eps(3);
            
            res.element(i).SX(j)  = sigma(1);
            res.element(i).SY(j)  = sigma(2);
            res.element(i).SXY(j) = sigma(3);
            
        end
        
    end

    % extract and store pure nodal results also
    xdof = 1:2:(2*nn-1);
    ydof = 2:2:2*nn;
    
    % displacements
    res.nodal.D   = D;
    res.nodal.UX  = D(xdof);
    res.nodal.UY  = D(ydof);
    res.nodal.U   = sqrt(D(xdof).^2 + D(ydof).^2);
    
    % constitutive matrix
    res.Em  = Em;

end


% Meshing
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
% [nodes] = node table: 1=x-coordinate, 2=y-coordinate, (row no. = node no.)
% [elems] = element table: 1=node i, 2=node j, ... (row no. = elem no.) 

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





%%%%%%%%%%%%%%%%%% plotting functions %%%%%%%%%%%%%%%%%%%%%
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
    nn = size(nodes,1);
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

function plot_results(elems,nodes,res,plot_mode)
    
    % add displacements to node coordinates to plot deformed shape
    UX = res.nodal.UX;
    UY = res.nodal.UY;
    nodes = nodes + [UX UY]; 
    
    % make vector of field quantity
    
       
    switch plot_mode
        case 'mesh' % just the deformed mesh
            patch('faces',elems,'vertices',nodes,'facecolor','none');
        
        case {'UX','UY','U'} % nodal results
            phi = eval(['res.nodal.' plot_mode]);
            patch('faces',elems,'vertices',nodes,'facecolor','interp','FaceVertexCData',phi);
        
        case {'SX','SY','SXY','EX','EY','EXY'} % element-based results
            for i=1:size(elems,1)
                phi = eval(['res.element(i).' plot_mode]);
                eface = 1:4;
                n = elems(i,1:4);
                enodes = nodes(n,:);
                patch('faces',eface,'vertices',enodes,'facecolor','interp','FaceVertexCData',phi');
            end
            
    end
    
    % enable colorbar
    if ~strcmp(plot_mode,'mesh')
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
