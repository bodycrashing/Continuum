close all

i = 9
    nc = i*2;
    no = nc;

%[nodes,elems] = Q4mesh_UnitCell(20,8,nc,no);
[nodes,elems,ndof] = MeshMaster('LSTiso',1,0.4,nc,no);

plot(nodes(:,1),nodes(:,2),'o')

%%
hold on
%test_x = reshape(nodes(:,1),nc+1,no+1);
%test_y = reshape(nodes(:,2),nc+1,no+1);

%elems = elems(:,[1,4,2,5,3,6]);

for i = 1:length(elems)

    x = nodes(elems(i,:),1);
    y = nodes(elems(i,:),2);
    
    x = [x;x(1)];
    y = [y;y(1)];
    
plot(x,y,'k')
plot(x',y','k')



end
axis equal
hold off


    % Model Plot
    figure()
    % (comment out what you don't need)
    init_plot(nodes); % always keep this one!
    plot_elements(elems,nodes)
    %plot_elem_numbers(elems,nodes)
    %plot_elem_cs(elems,nodes)
    plot_nodes(nodes)
    plot_node_numbers(nodes)
    %plot_dof_numbers(nodes)
    %plot_BCs(nodes,BC)
    %plot_loads(nodes,P)
    %plot_dist_load(elems,nodes,poly,W)
    
%%


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
            text(nodes(i,1)+0.01,nodes(i,2)-0.01,['(' num2str(i*2-1) ',' num2str(i*2) ')'],'color','b','fontsize',7);
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
        
        if size(elems,2)==6
            elems = elems(:,[1,4,2,5,3,6]);
        end
        
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


