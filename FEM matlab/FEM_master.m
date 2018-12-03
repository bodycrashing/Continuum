function [res,E_unitCell,nu_unitCell,SigmaVM_max,n_elems,eta] = FEM_master(f,meshType,nc,no,plotTF,mesher)
% f = porousity
% meshType = Can be 'CSTiso', 'LSTiso', 'Q4iso'
% nc = number of elements allong hole curve
% no = number of elements allong symetry boundary
% plotTF = true/false specifies if the function should create plots
% Optional :
% mesher = if mesher is different from meshType. This can be:
% 'FuldPladeLST' or 'FuldPladeQ4'.

% If no special mesher is defined:
if nargin < 6
    mesher = meshType;
end


%% %%%%%%%%%%%%%%%%%% define model %%%%%%%%%%%%%%%%%%%%%%%%%
W  = 1;        % Unit zell widht and height [mm]
rad  = 2*W*sqrt(f/(pi));        % Radius of hole [mm]
h  = 1;         % model thickness
%rad = 8/10*W;
nGP = 2;        % no. Gauss points

% Material properties
E  = 2230.6;    % Youngs modulus
nu = 0.355;       % Poissons ratio

D_forced = 1/W; %Forced displacement [mm]

%% %%%%%%%% Generate Mesh, Constraints and Loads %%%%%%%%%%%%%%%%%%%%%%%

% generate mesh (node and element tables)
[nodes,elems,ndof] = MeshMaster(mesher,W,rad,nc,no);

% set boundary conditions
BC = false(ndof,1);
MC = false(ndof,1);
% fixed dofs (at x = 0 and y = 0), all others are free
BC(find(nodes(:,1) == 0)*2-1) = 1; % Left side nodes
BC(find(nodes(:,2) == 0)*2) = 1; % Bottom nodes
%BC(find(nodes(:,1) == W)*2-1) = 1; % Right side nodes
BC(find(nodes(:,2) == 1)*2) = 1; % Top nodes
% Forced displacements
Dc = zeros(ndof,1);
Dc(find(nodes(:,2) == 1)*2) = D_forced;

% Multivariable constraints
MC(find(nodes(:,1) == 1)*2-1) = 1; % =1 if on
%MC(find(nodes(:,1) == -1.5)*2-1) = 1; % =1 if on
C = multi_constraints(elems,nodes,MC);

P = zeros(ndof,1);
R = P;


%% %%%%%%%%%%%%%%%%%%%% setup & solve %%%%%%%%%%%%%%%%%%%%%%%
% constitutive matrix for plane stress, Cook eq. (3.1-3)
Em = [ 1 nu 0;
    nu  1 0;
    0  0 (1-nu)/2]*E/(1-nu^2);

% Assembly of global stiffness matrix
K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP,meshType);

% solve for displacements and reaction
[D,R] = solve(K,R,BC,Dc,C);

% Calculating stresses in the elements
res = post_process(elems,nodes,D,Em,meshType);

eta = estimate_error(res,elems,nodes,meshType);

F_res = sum(R(find(nodes(:,2) == 1)*2));

D_y = (D(find(nodes(:,1) == 1,1)*2-1));

E_unitCell = (F_res)/(D_forced);

nu_unitCell = -D_y/D_forced;


for i = 1:length(res.element)
    Sig(:,i) = res.element(i).Se;
end
SigmaVM_max = max(Sig(:));
n_elems = length(elems);


%% %%%%%%%%%%%%%%%%%%%% plot model %%%%%%%%%%%%%%%%%%%%%%%%%%

if plotTF  
    % Model Plot
    % (comment out what you don't need)
    init_plot(nodes); % always keep this one!
    plot_elements(elems,nodes)
    %plot_elem_numbers(elems,nodes)
    %plot_elem_cs(elems,nodes)
    plot_nodes(nodes)
    plot_node_numbers(nodes)
    %plot_dof_numbers(nodes)
    plot_BCs(nodes,BC)
    plot_loads(nodes,P)
    %plot_dist_load(elems,nodes,poly,W)
    plot_reaction_forces(nodes,D,R,BC)
    
    figure()
    % Results Plot
    plot_results(elems,nodes,res,'Seavg')
    
    figure()
    plot_displacement(elems,nodes,D,'mesh') %last argument: 'UX', 'UY', 'U' or 'mesh'
    
end



%% %%%%%%%%%%%%%%%% FE functions %%%%%%%%%%%%%%%%%%%%%%%%%%%

    function k = k_elem(x_vec,y_vec,Em,h,nGP,meshType)
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
        k = zeros(length(x_vec)*2);
        % If statement determines if triangle or square element
        if mod(length(x_vec),3) ~= 0 % Rectangular element
            % Get Gauss points and weights
            [GP,GPW] = Gauss(nGP);
            
            for i = 1:nGP
                for j = 1:nGP     
                    % xi and eta are set to be the gaus points
                    xi = GP(i);
                    eta = GP(j);
                    % Evaluating the B matrix at xi and eta
                    [B,J_det] = BmatMaster(meshType,x_vec,y_vec,xi,eta);
                    
                    k = k + B'*Em*B*h*abs(J_det)*GPW(i)*GPW(j);
                end
            end    
        else % Triangular element
            Int_points = [2/3, 1/6
                          1/6, 1/6
                          1/6, 2/3];
            
            GPW = [1/3, 1/3, 1/3];  
            for i = 1:3
                [B,J_det] = BmatMaster(meshType,x_vec,y_vec,Int_points(i,1),Int_points(i,2));
                
                k = k + B'*Em*B*h*(1/2)*abs(J_det)*GPW(i);
            end
        end
    end






    function K = global_stiffness_matrix(elems,nodes,Em,h,ndof,nGP,meshType)
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

        K = zeros(ndof,ndof);    
        for i = 1:size(elems,1)
            nod = elems(i,:);
            
            x_vec = nodes(nod,1);
            y_vec = nodes(nod,2);
            
            dof_x = nod*2-1;
            dof_y = nod*2;
            dof = reshape([dof_x;dof_y],1,[]);
            
            k = k_elem(x_vec,y_vec,Em,h,nGP,meshType);
            
            K(dof,dof) = K(dof,dof) + k;
            
        end 
    end



    function [GP,GPW] = Gauss(nGP)
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
        GPW = w_l(nGP,1:nGP);
    end


    function [D,R] = solve(K,R,BC,Dconstrains,C)
        % solve FE problem, Cook eq. (2.7-3, 2.7-2b)
        % D = vector of displacements (both free and fixed)
        % R = vector of nodal forces (applied and reactions)
        
        O = zeros(size(C,1));
        Q = zeros(size(C,1),1);
        
        BC2 = cat(1,BC,false(size(C,1),1));
        
        ConstraintMat = [K C'
                         C O];
        
        % Splitting global stiffness matrix
        K11 = ConstraintMat(~BC2,~BC2);
        K12 = ConstraintMat(~BC2,BC2);
        K21 = ConstraintMat(BC2,~BC2);
        K22 = ConstraintMat(BC2,BC2);
        
        % Defining the constrained deformation and forces
        Dc = Dconstrains(BC);
        Rc = R(~BC);
        
        Rc = cat(1,Rc,Q);
        
        % Calculating free deformations and forces
        Dx = K11\(Rc-K12*Dc);
        Rx = K21*Dx + K22*Dc;
        
        % Filling D and R vector
        D = zeros(length(BC),1);
        R = zeros(length(BC),1);
        
        D(~BC) = Dx(1:end-length(Q));
        D(BC) = Dc;
        
        R(~BC) = Rc(1:end-length(Q));
        R(BC) = Rx;
        
        
    end


    function C = multi_constraints(elems,nodes,MC)
        test = [1 0 -1];
        ndof = length(nodes)*2;
        C = zeros(sum(MC)-1,ndof);
        
        for i = 1:sum(MC)-1
            C(i,(i*2-1:i*2+1)) = test;
        end
        
    end


    function res = post_process(elems,nodes,D,Em,meshType)
        
        % numbers of elements and nodes
        ne = size(elems,1);
        nn = size(nodes,1);
        
        % loop over elements: calculate nodal results for each element (unaveraged)
        for i = 1:ne
            
            % nodes and dof's for element i
            n = elems(i,:);
            
            dof_x = n*2-1;
            dof_y = n*2;
            dof = reshape([dof_x;dof_y],1,[])';
            
            d = D(dof); % element displacement vector
            
            % global node coordinates
            x = nodes(n,1);
            y = nodes(n,2);
            
            % local node coordinates
            if strcmp(meshType,'Q4iso')
                c = [-1 -1;
                    1 -1;
                    1  1;
                    -1  1];
                scale = sqrt(3);
            elseif strcmp(meshType,'CSTiso') || strcmp(meshType,'LSTiso')
                
                c = [0,0;
                    1/2,0;
                    1,0;
                    1/2,1/2;
                    0,1;
                    0,1/2];
                scale = 1;
            end
            
            % loop over element nodes
            for j = 1:size(elems,2)
                
                % evaluate B matrix at node j
                [B,~,~] = BmatMaster(meshType,x,y,scale*c(j,1),scale*c(j,2));
                
                % calculate strain and stress the simple way:
                %eps   = B*d;
                %sigma = Em*eps;
                epsGP(:,j)   = B*d;
                sigmaGP(:,j) = Em*epsGP(:,j);
                
            end
            
            for j = 1:size(elems,2)
                [~,~,N] = BmatMaster(meshType,x,y,scale*c(j,1),scale*c(j,2));
                
                eps = N*epsGP';
                sigma = N*sigmaGP';
                
                Sx  = sigma(1);
                Sy  = sigma(2);
                Txy = sigma(3);
                
                % von Mises
                Se  = sqrt( Sx^2 + Sy^2 -Sx*Sy + 3*Txy^2 );
                
                % principal stresses
                S1  = (Sx+Sy)/2 + sqrt( ((Sx-Sy)/2)^2 + Txy^2);
                S2  = (Sx+Sy)/2 - sqrt( ((Sx-Sy)/2)^2 + Txy^2);
                
                % store results per element
                res.element(i).EX(j)  = eps(1);
                res.element(i).EY(j)  = eps(2);
                res.element(i).EXY(j) = eps(3);
                
                res.element(i).SX(j)  = Sx;
                res.element(i).SY(j)  = Sy;
                res.element(i).SXY(j) = Txy;
                
                res.element(i).S1(j)  = S1;
                res.element(i).S2(j)  = S2;
                res.element(i).Se(j)  = Se;
                
            end
            
        end
        
        % Stress averaging
        for j = 1:max(elems(:)) % Looping over number of nodes
            [elem_nr,node_pos] = find(elems==j); % Finding element nr and pos for each node
            node_Se = [];
            node_SX = [];
            node_SY = [];
            node_SXY = [];
            for i = 1:length(elem_nr)
                % Finding stresses for the node
                node_Se(i) =  res.element(elem_nr(i)).Se(node_pos(i));
                node_SX = res.element(elem_nr(i)).SX(node_pos(i));
                node_SY = res.element(elem_nr(i)).SY(node_pos(i));
                node_SXY = res.element(elem_nr(i)).SXY(node_pos(i));
            end
            % Averaging the stresses
            node_Se_avg = mean(node_Se);
            node_SX_avg = mean(node_SX);
            node_SY_avg = mean(node_SY);
            node_SXY_avg = mean(node_SXY);
            
            % Inserting the stresses back in the global stress matrix
            for i = 1:length(elem_nr)
                res.element(elem_nr(i)).Seavg(node_pos(i)) = node_Se_avg ;
                res.element(elem_nr(i)).SXavg(node_pos(i)) = node_SX_avg ;
                res.element(elem_nr(i)).SYavg(node_pos(i)) = node_SY_avg ;
                res.element(elem_nr(i)).SXYavg(node_pos(i)) = node_SXY_avg ;
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



    function R = element_loads(poly,elems,nodes,h,W,nGP)
        % !!!!!!!!!!! THIS FUNCTION ONLY WORK ON Q4 !!!!!!!!!!
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
                    N1 = 1/(4*1*1)*(1-xi)*(1-eta);
                    N2 = 1/(4*1*1)*(1+xi)*(1-eta);
                    N3 = 1/(4*1*1)*(1+xi)*(1+eta);
                    N4 = 1/(4*1*1)*(1-xi)*(1+eta);
                    
                    N = kron([N1 N2 N3 N4],eye(2));
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



    function eta = estimate_error(res,elems,nodes,meshType)
        % Following cook 9.10
        
         if false%mod(length(x_vec),3) ~= 0
             [Int_points,GPW] = Gauss(2);
         else
            Int_points = [2/3, 1/6
                          1/6, 1/6
                          1/6, 2/3];
              GPW = [1/3, 1/3, 1/3];
         end
              
              U = 0;
              e = 0;
        for i = 1:length(elems)
            nod = elems(i,:);
            
            x_vec = nodes(nod,1);
            y_vec = nodes(nod,2);
            Ui = 0;
            ei = 0;
            % Integrating using Gauss Quadrature
            for j = 1:3
            [~,j_det,N] = BmatMaster(meshType,x_vec,y_vec,Int_points(j,1),Int_points(j,2));
            
            sigma = (N* [res.element(i).SX' , res.element(i).SY' , res.element(i).SXY'])';
            
            sigma_avg = (N * [res.element(i).SXavg' , res.element(i).SYavg' , res.element(i).SXYavg'])';
            
            Ui = Ui + sigma'*sigma*1/2*abs(j_det)*GPW(j);
            
            ei = ei + (sigma_avg-sigma)' * (sigma_avg-sigma) *1/2*abs(j_det)*GPW(j);
            end
        U = U + Ui;
        e = e + ei;
        
        eta  = e / (U + e);
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
                        nvec' + 0.3*[-0.5  0.25];
                        nvec' + 0.3*[-0.5 -0.25]];
                    patch('faces',[1 2 3],'vertices',v_verts,'edgecolor','k','facecolor','none','tag','Boundary conditions');
                else % y-dof
                    h_verts = [nvec';
                        nvec' + 0.3*[0.25 -0.5];
                        nvec' + 0.3*[-0.25 -0.5]];
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
            text(nodes(i,1)+0.01,nodes(i,2)+0.01,num2str(i),'color','r');
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
            
            %elems = elems(:,[1,4,2,5,4,6]);
            
            % nodes on element
            n = elems(i,:);
            
            % element retangles
            verts = nodes(n,:);
            
            if size(elems,2)==6
                verts = verts([1,4,2,5,3,6],:);
            end
            
            patch('faces',1:length(n),'vertices',verts,'facecolor','c','edgecolor','k');
            
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
            n = elems(i,:);
            
            % element CS & number
            mid = [0 0]';
            for j = 1:length(n)
                vec = [nodes(n(j),1) nodes(n(j),2)]';  % vectors to nodes
                mid = mid + vec;
            end
            mid = mid/length(n);
            
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
        % i is the element number
        % j is the node number in each element
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
                sigma_vm = real(sigma_vm);
            end
            if size(elems,2)==6
                enodes = enodes([1,4,2,5,3,6],:);
                sigma_vm = sigma_vm([1,4,2,5,3,6],:);
            end
            
            h(i) = patch('faces',eface,'vertices',enodes,...
                'facecolor','interp','FaceVertexCData',sigma_vm(:,i));
            
            %h(i) = patch('faces',eface,'vertices',enodes,...
            %     'facecolor','interp','FaceVertexCData',sigma_val(1,:,i)');
            
            
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
        
        colormap(jet(10))
        
    end




    function plot_results(elems,nodes,res,plot_mode)
        
        % add displacements to node coordinates to plot deformed shape
        UX = res.nodal.UX;
        UY = res.nodal.UY;
        %nodes = nodes + [UX UY];
        
        switch plot_mode
            case 'mesh' % just the deformed mesh
                patch('faces',elems,'vertices',nodes,'facecolor','none');
                
            case {'UX','UY','U'} % nodal results
                phi = eval(['res.nodal.' plot_mode]);
                patch('faces',elems,'vertices',nodes,'facecolor','interp','FaceVertexCData',phi);
                
            case {'SX','SY','SXY','EX','EY','EXY' 'S1' 'S2' 'Seavg'} % element-based results
                for i = 1:size(elems,1)
                    phi = eval(['res.element(i).' plot_mode]);
                    eface = 1:size(elems,2);
                    n = elems(i,:);
                    enodes = nodes(n,:);
                    
                    if size(elems,2)==6
                        enodes = enodes([1,4,2,5,3,6],:);
                        phi = phi([1,4,2,5,3,6]);
                    end
                    
                
                    patch('faces',eface,'vertices',enodes,'facecolor','interp','FaceVertexCData',phi');
                end
                
        end
        
        % enable colorbar
        if ~strcmp(plot_mode,'mesh')
            cb = colorbar;
            colormap(jet(9));
            cb.Title.String = plot_mode;
        end
        
    end

end