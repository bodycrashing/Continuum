function [nodes,elems,ndof] = T3mesh(W,R,no,nc)
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

r = linspace(R/W,1,no+1); % Vi bygger her dobbelt så mange noder som i Q4 modellen!
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
