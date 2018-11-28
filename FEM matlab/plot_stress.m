function plot_stress(results,stress_type,average,plot_full)
% This function generates a stress plot of the entire RVE
%
% INPUTS:
%       results = structure containing FE problem results
%   stress_type = string detailing which stress to plot ('x','y','xy','von')
%       average = logical operator telling whether to plot average or
%                 unaveraged nodal stress
%     plot_full = boolean, 0 only plots a quarter of the RVE, 1 plots all
    
    % initialize figure
    createFigure(0.4,0.7);
    axis equal
    hold on
    
    % extract element and node table
    elems = results.elems;
    nodes = results.nodes;
    
    % no. elements
    ne = size(elems,1);
    
    % are we averaging the stresses?
    if average
        
        % get average stress on all nodes
        aveStress = averageStress(results);
        
        % initialize zero matrix
        stress = zeros(ne,size(elems,2),3);
        for i = 1:ne
            for j = 1:size(elems,2)
                
                % assign node stress for each element
                stress(i,j,:) = aveStress(elems(i,j),:);
            
            end
        end
        
        switch stress_type
            case 'x'
                stress = stress(:,:,1);
            case 'y'
                stress = stress(:,:,2);
            case 'xy'
                stress = stress(:,:,3);
            case 'von'
                sig_x  = stress(:,:,1);
                sig_y  = stress(:,:,2);
                sig_yx = stress(:,:,3);
                
                % Von Mises stress for general plane stress
                stress = sqrt((sig_x-sig_y).^2 + 3*sig_yx.^2);
        end
        
    else
        % extract the wanted stress data
        switch stress_type
            case 'x'
                stress = results.stress(:,:,1);
            case 'y'
                stress = results.stress(:,:,2);
            case 'xy'
                stress = results.stress(:,:,3);
            case 'von'
                sig_x  = results.stress(:,:,1);
                sig_y  = results.stress(:,:,2);
                sig_yx = results.stress(:,:,3);
                
                % Von Mises stress for general plane stress
                stress = sqrt(sig_x.^2 - sig_x.*sig_y + sig_y.^2 + 3*sig_yx.^2);
        end
    end
    
    % loop over every element and plot a patch for each
    if plot_full
        for i = -1:2:1
            for j = -1:2:1
                Nodes = [nodes(:,1)*i nodes(:,2)*j];
                for k = 1:ne
                    if size(elems,2) == 9 % Q9 element (because it has a midnode that we don't want to plot)
                        patch('faces',1:8,'vertices',Nodes(elems(k,1:8),:),'EdgeColor','k','facecolor','interp','FaceVertexCData',stress(k,1:8)');
                    else % every other element
                        patch('faces',1:size(elems,2),'vertices',Nodes(elems(k,:),:),'EdgeColor','k','facecolor','interp','FaceVertexCData',stress(k,:)');
                    end
                end
            end
        end
    else
        for i = 1:ne
            if size(elems,2) == 9 % Q9 element (because it has a midnode that we don't want to plot)
                patch('faces',1:8,'vertices',nodes(elems(i,1:8),:),'EdgeColor','k','facecolor','interp','FaceVertexCData',stress(i,1:8)');
            else % every other element
                patch('faces',1:size(elems,2),'vertices',nodes(elems(i,:),:),'EdgeColor','k','facecolor','interp','FaceVertexCData',stress(i,:)');
            end
        end
    end
    
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    cb.Title.String = '[MPa]';
    colormap(jet(9));
    
    xlabel('x [mm]')
    ylabel('y [mm]')
    title(['$\sigma_{',stress_type,'}$'],'interpreter','latex')
end


function S = averageStress(results)
% This function calculates the average stress for all nodes in "results".
%
% INPUT:
% results = structure containing FE problem results
%
% OUTPUT:
%       S = average Von Mises stress for all nodes

    % extract element and node table
    elems = results.elems;
    nodes = results.nodes;
    
    % no. elements
    ne = size(elems,1);
    
    % no. nodes
    nn = size(nodes,1);
    
    % extract stresses
    sigma = results.stress;
    
    % initialize average stress table
    S = zeros(nn,3);
    n = zeros(nn,1);
    
    for i = 1:ne
        for j = 1:size(elems,2)
            
            % get current node
            cn = elems(i,j);
            
            % reshape 1x1x3 to a 1x3 vector
            stresses = reshape(sigma(i,j,:),[1 3]);
            
            % increment nodal stress
            S(cn,:) = S(cn,:) + stresses;
            
            % increment counter
            n(cn) = n(cn) + 1;
            
        end
    end
    
    % divide nodal stresses with the stress counter
    S = S./n;

end