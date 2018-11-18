function [nodes,elems,ndof] = MeshMaster(type,W,R,nc,no)
% type = string describing element type.
% possible values: 'Q4' 'Q8' 'CST' 'LST'
%
% nc = number of elements on hole curve
% no = number of elements in direction normal to hole curve
% R = Radius of inner holea
% W = Width of block

% Include mesh functions under suitable cases.
switch type
    case 'Q4iso'
        s_increase = 1.3; % Speed which the superelipse becomse square
        % no. of dof's
        ndof = 2*(nc+1)*(no+1);
        
        % Vector containing radiuses
        r = linspace(R/W,1,no+1);
        % Vector containing node positions on boundary
        b = linspace(0,1,nc/2+1);
        
        % Creating boundary nodes
        for i = 1:((nc+1)/2)+1
            if i == 1
                nodes(i,1) = 1;
                nodes(i,2) = 0;
            elseif i == ((nc+1)/2)+1
                nodes(i,1) = 1;
                nodes(i,2) = 1;
            else
                nodes(i,1) = 1;
                nodes(i,2) = b(i);
            end
        end
        % Mirroring the created nodes arround diagonal
        nodes = cat(1,nodes,rot90(nodes(1:end-1,:),2));
        
        n = 0; % Node number
        
        % Creating nodes following superelipse
        for k = no:-1:1
            nodes2=[]; % Temporary node container
            for i = 0:(nc/2)
                n = n+1; % Increase node number
                
                n_elep = 2+(k-1)/s_increase; % Order of Superelipse
                
                angle = pi/4 .* (i/(nc/2)).^(n_elep/2); % Angle to evaluate superelipse at
                
                % Generating node x and y coordinate by evaluating the superelipse
                % at angele
                nodes2(i+1,1) = r(k) * cos(angle) ^(2/n_elep);
                nodes2(i+1,2) = r(k) * sin(angle) ^(2/n_elep);
                
            end
            % Mirroring nodal points arround diagonal and concattinating with
            % global node vector
            nodes2 = cat(1,nodes2,rot90(nodes2(1:end-1,:),2));
            nodes = cat(1,nodes,nodes2);
        end
        
        % Generating book keeping matrix
        e_no = 0;
        for j=1:no
            for i = 1:nc
                e_no = e_no+1;
                
                n1 = i + (nc+1)*j - (nc+1);
                n2 = n1 + nc+1;
                n3 = n2 + 1;
                n4 = n1 + 1;
                
                elems(e_no,1:4) = [n1 n2 n3 n4];
            end
        end
        
        nodes = nodes * W;
        elems = fliplr(elems);
        elems =  elems(:,[3 4 1 2]);
        
        
    case 'Q8'
        disp('No Q8 function yet')
        
    case 'CSTiso'
        %[nodes,elems,ndof] = ImportMesh('CST2_nodes.txt','CST2_elems.txt');
        if true
            % Square panel mesh generation with 3-noded triangle elements (T3/CST).
            %
            % INPUTS:
            % W  = width of model
            % H  = height of model
            % nw = no. of elements in x-direction
            % nh = no. of elements in y-direction
            %
            % OUTPUTS:
            % nodes = node table: 1=x-coordinate, 2=y-coordinate, (row no. = node no.)
            % elems = element table: 1=node i, 2=node j, ... (row no. = elem no.)
            
            s_increase = 1.3; % Speed which the superelipse becomse square
            
            % no. of dof's
            ndof = 2*(no+1)*(nc+1);
            
            % no. of nodes
            nn = ndof/2;
            
            % no. elements
            ne = 2*no*nc;
            
            % initialize array
            elems = zeros(ne,3);
            
            r = linspace(R/W,1,no+1); % Vi bygger her dobbelt s? mange noder som i Q4 modellen!
            % Vector containing node positions on boundary
            b = linspace(0,1,ceil((nc+1)/2));
            
            % Creating boundary nodes
            for i = 1:length(b)
                if i == 1
                    nodes(i,1) = 1;
                    nodes(i,2) = 0;
                elseif i == length(b)
                    nodes(i,1) = 1;
                    nodes(i,2) = 1;
                else
                    nodes(i,1) = 1;
                    nodes(i,2) = b(i);
                end
            end
            % Mirroring the created nodes arround diagonal
            nodes = cat(1,nodes,rot90(nodes(1:end-1,:),2));
            
            
            % Creating nodes following superelipse
            for k = no:-1:1
                nodes2=[]; % Temporary node container
                for i = 1:length(b)
                    
                    n_elep = 2+(k-1)/s_increase; % Order of Superelipse
                    
                    angle = pi/4 .* ((i-1)/(nc/2)).^(n_elep/2); % Angle to evaluate superelipse at
                    
                    % Generating node x and y coordinate by evaluating the superelipse
                    % at angele
                    nodes2(i,1) = r(k) * cos(angle) ^(2/n_elep);
                    nodes2(i,2) = r(k) * sin(angle) ^(2/n_elep);
                    
                end
                % Mirroring nodal points arround diagonal and concattinating with
                % global node vector
                nodes2 = cat(1,nodes2,rot90(nodes2(1:end-1,:),2));
                nodes = cat(1,nodes,nodes2);
            end
            
            % setup element table
            e = 0;
            for j=1:no
                for i = 1:nc
                    e = e+1;
                    
                    n1 = (nc+1)*j + i -(nc+1);
                    n2 = (nc+1) + n1;
                    n3 = n2 + 1;
                    n4 = n1 + 1;
                    
                    if i <= ceil(nc/2)
                        elems(e,1:3) = [n1 n2 n4];
                        e = e+1;
                        elems(e,1:3) = [n3 n4 n2];
                    else
                        elems(e,1:3) = [n1 n2 n3];
                        e = e+1;
                        elems(e,1:3) = [n3 n4 n1];
                    end
                end
            end
            nodes = nodes * W;
        end
        
    case 'LSTiso'
        
        % The function generates a mesh for a unit cell with width W and hole
        % radius R
        %
        % nc = number of elements on hole curve
        % no = number of elements in direction normal to hole curve
        % R = Radius of inner holea
        % W = Width of block
        
        s_increase = 1.2; % Speed which the superelipse becomse square
        
        % no. of dof's
        nc = 2*nc; % Number of nodes along the cirumfrance of the curve
        no = 2*no; % Number of nodes
        
        % Vector containing radiuses
        
        r = linspace(R/W,1,no+1); % Vi bygger her dobbelt s? mange noder som i Q4 modellen!
        % Vector containing node positions on boundary
        b = linspace(0,1,ceil((nc+1)/2));
        
        % Creating boundary nodes
        for i = 1:length(b)
            if i == 1
                nodes(i,1) = 1;
                nodes(i,2) = 0;
            elseif i == length(b)
                nodes(i,1) = 1;
                nodes(i,2) = 1;
            else
                nodes(i,1) = 1;
                nodes(i,2) = b(i);
            end
        end
        % Mirroring the created nodes arround diagonal
        nodes = cat(1,nodes,rot90(nodes(1:end-1,:),2));
        
        
        % Creating nodes following superelipse
        for k = no:-1:1
            nodes2=[]; % Temporary node container
            for i = 1:length(b)
                
                n_elep = 2+(k-1)/s_increase; % Order of Superelipse
                
                angle = pi/4 .* ((i-1)/(nc/2)).^(n_elep/2); % Angle to evaluate superelipse at
                
                % Generating node x and y coordinate by evaluating the superelipse
                % at angele
                nodes2(i,1) = r(k) * cos(angle) ^(2/n_elep);
                nodes2(i,2) = r(k) * sin(angle) ^(2/n_elep);
                
            end
            % Mirroring nodal points arround diagonal and concattinating with
            % global node vector
            nodes2 = cat(1,nodes2,rot90(nodes2(1:end-1,:),2));
            nodes = cat(1,nodes,nodes2);
        end
        
        % Generating book keeping matrix
        e_no = 0;
        for j=1:no/2
            for i = 1:nc/2
                
                e_no = e_no + 1;
                
                q = 2*i-1;
                
                n1 = 2*(nc+1)*j + q -2*(nc+1);
                n2 = (nc+1) + n1;
                n3 = (nc+1) + n2;
                n4 = n1 + 1;
                n5 = n2 + 1;
                n6 = n3 + 1;
                n7 = n4 + 1;
                n8 = n5 + 1;
                n9 = n6 + 1;
                
                if i <= ceil(nc/4)
                    elems(e_no,1:6) = [n1 n2 n3 n5 n7 n4];
                    e_no = e_no + 1;
                    elems(e_no,1:6) = [n9 n8 n7 n5 n3 n6];
                else
                    elems(e_no,1:6) = [n1 n2 n3 n6 n9 n5];
                    e_no = e_no + 1;
                    elems(e_no,1:6) = [n9 n8 n7 n4 n1 n5];
                end
            end
        end
        
        elems = fliplr(elems);
        elems = elems(:,[6,2,4,1,3,5]);
        nodes = nodes * W;
        ndof = 2*max(elems(:));
        
        %[nodes,elems,ndof] = ImportMesh('LST_nodes.txt','LST_elems.txt');
end


    function [nodes,elems,ndof] = ImportMesh(file_nodes,file_elem)
        
        nodes_t = readtable(file_nodes);
        elems_t = readtable(file_elem);
        
        nodes = table2array(nodes_t(:,2:3));
        nodes = nodes*1000;
        elems = table2array(elems_t(:,3:size(elems_t,2)-1));
        
        ndof = length(nodes)*2;
        
        
    end
end