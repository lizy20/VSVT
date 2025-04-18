function [ nodecoor_list, nodecoor_cell, ele_cell ] = getNodeEle( vert, tria, tnum, ele_order )
% get node coordinares, and elements
% input:
%   vert(k,1:2) = [x_coordinate, y_coordinate] of k-th node 
%   tria(m,1:3) = [node_numbering_of_3_nodes] of m-th element
%   tnum(m,1) = n; means the m-th element is belong to phase n
% output:
%   nodecoor_cell{i} and ele_cell{i} represent phase i
%   nodecoor_cell{i}(j,1:3) = [node_numbering, x_coordinate, y_coordinate]
%   ele_cell{i}(j,1:4/7) = [element numbering, node_numbering_of_all_nodes]
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Oct 2020
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   

    if ele_order == 1
        % do nothing
    elseif ele_order == 2
        % insert midpoint in each edges
        % linear element to quadratic element
        [ vert, tria ] = insertNode( vert, tria );
        % size of tria will change after insertNode
        % tria(m,1:6) = [node_numbering_of_6_nodes] of m-th element
    else
        error('ele_order should be 1 or 2')
    end
    
    nodecoor_list = zeros( size(vert,1), 3 );
    nodecoor_list( :, 1 ) = 1: size(vert,1);
    nodecoor_list( :, 2:3 ) = vert;
    
    % extract the mesh of each phase from meshes
    [ nodecoor_cell, ele_cell ] = whole2phase( vert,tria,tnum );
    
end

function [ nodecoor_cell, ele_cell ] = whole2phase( vert, tria, tnum )
% whole2phase: extract the mesh of each phase from meshes
% input:
%   vert(k,1:2) = [x_coordinate, y_coordinate] of k-th node 
%   tria(m,:) = [node_numbering_of_all_nodes] of m-th element (3 or 6 nodes)
%   tnum(m,1) = n; means the m-th element is belong to phase n
% output:
%   nodecoor_cell{i} and ele_cell{i} represent phase i
%   nodecoor_cell{i}(j,1:3) = [node_numbering, x_coordinate, y_coordinate]
%   ele_cell{i}(j,1:4/7) = [element numbering, node_numbering_of_all_nodes]
%

    phase_code_vec = unique( tnum );
    num_phase = length( phase_code_vec );
    nodecoor_cell = cell( 1, num_phase );
    ele_cell = cell( 1, num_phase );

    for i = 1: num_phase
        % get triangular mesh (node numbering of 3 nodes) of phase_i from tria
        tria_phase_i = tria( tnum==phase_code_vec(i), : );
        
        % get unique node numbering of phase_i
        nodes_phase_i = unique( tria_phase_i(:) );
        
        % extract node coordinates from vert according to node numbering
        nodecoor_phase_i = zeros( length(nodes_phase_i), 3 );
        nodecoor_phase_i(:,1) = nodes_phase_i;  % node numbering
        for j = 1: size(nodecoor_phase_i,1)
              nodecoor_phase_i(j,2:3) = vert( nodecoor_phase_i(j,1), : );
        end
        
        % get element numbering and node numbering of phase_i
        ele_phase_i = zeros( size(tria_phase_i,1), 1+size(tria,2) );
        % element numbering
        ele_phase_i(:,1) = find( tnum==phase_code_vec(i) );
        % node numbering of 3 nodes or 6 nodes
        ele_phase_i(:,2:end) = tria_phase_i;
        
        nodecoor_cell{i} = nodecoor_phase_i;
        ele_cell{i} = ele_phase_i;
    end
end


function [ vert, tria ] = insertNode( vert, tria )
% insertNode: insert a midpoint in each edges,
%            linear element -> quadratic element (like CPS6M in abaqus)

    num_node = size(vert,1);
    num_ele = size(tria,1);
    
    % append first index to the end
    tria_app = [ tria, tria(:,1) ];

    for j = 1: num_ele
        for k=1:3
            % index of node numbering
            idx1 = tria_app(j,k);
            idx2 = tria_app(j,k+1);
            % get coor
            coor1 = vert( idx1, 1:2 );
            coor2 = vert( idx2, 1:2 );
            % get middle point
            coor3 = (coor1 + coor2)/2;
            
            num_node = num_node + 1;
            % add coordinates of new vertex into vert
            vert( num_node, 1:2 ) = [ coor3(1) coor3(2) ];
            % add node numbering of new vertex into tria
            tria( j, k+3 ) = num_node;
        end
    end
end
