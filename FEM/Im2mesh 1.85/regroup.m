function [ nodeU, edgeU, part ] = regroup( node_cell, edge_cell )
% group node_cell, edge_cell into nodeU, edgeU, part

    % ---------------------------------------------------------------------
    % group node_cell, edge_cell into node, edge, part
    
    % node_cell and edge_cell have the same size
    num_phase = length( node_cell );
    
    % initialize edge and node
    % pre-allocated memory
    total_size = 0;
    for i = 1: num_phase
        total_size = total_size + size( node_cell{i}, 1 );
    end
    edge = zeros( total_size, 3 );
    node = zeros( total_size, 2 );
    
    accumu_size = 0;
    for i = 1: num_phase
        edge_cell{i}(:,3) = i;  % edge_cell{i}(:,3) as marker, show phase
        if i > 1
            edge_cell{i}(:,1:2) = edge_cell{i}(:,1:2) + accumu_size;
        end
        
        idx = (accumu_size + 1) : (accumu_size + size(node_cell{i}, 1));
        edge( idx, : ) = edge_cell{i};
        node( idx, : ) = node_cell{i};
        
        accumu_size = accumu_size + size( node_cell{i}, 1 );
    end

    % above code equal to:
    %     edg1(:,3) = 1;
    %     edg2(:,3) = 2;
    %     edg3(:,3) = 3;
    %     edg2(:,1:2) = edg2(:,1:2) + size(nod1,1);
    %     edg3(:,1:2) = edg3(:,1:2) + size(nod1,1) + size(nod2,1);
    %     edge = [edg1; edg2; edg3];
    %     node = [nod1; nod2; nod3];
	
    %-- the PART argument is a cell array that defines individu-
    %-- al polygonal "parts" of the overall geometry. Each elem-
    %-- ent PART{I} is a list of edge indexes, indicating which
    %-- edges make up the boundary of each region.

    part = cell( 1, num_phase );
    for i = 1: num_phase
        part{i} = find( edge(:,3) == i );
    end
        
    edge = edge( :, 1:2 );
    
    % ---------------------------------------------------------------------
    % Numbering of some nodes and edges aren't unique.
    % Here, we do re-numbering so they would be unique.
    
    % unique node
    [nodeU,~,icn] = unique(node,'rows');    % nodeU = icn(node)
    
    % replace edge using new index
    for i = 1: size(edge,1)
        edge(i,1) = icn( edge(i,1) );
        edge(i,2) = icn( edge(i,2) );
    end
    
    % unique edge
    edge=sort(edge,2);
    [edgeU,~,ice] = unique(edge,'rows');    % edgeU = ice(edge)

    % replace part using new index
    % part is edge index
    for i = 1: length(part)
       for j = 1: length(part{i})
           part{i}(j) = ice( part{i}(j) );
       end
    end
    % ---------------------------------------------------------------------
end

