
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




