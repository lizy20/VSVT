function showInterfNode( interf_node_part,nodecoor_part )
% showInterfNode: use to check whether interfacial nodes have been found
%                successfullly
% example
% [ interfBC_node_B, interfBC_node_C ] = findInterfNode( nodecoor_B, nodecoor_C );
% showInterfNode( interfBC_node_B, nodecoor_B );

    for i=1:length(interf_node_part)
        idx = interf_node_part(i);
        x = nodecoor_part( nodecoor_part(:,1)==idx, 2 );
        y = nodecoor_part( nodecoor_part(:,1)==idx, 3 );
        plot(x,y,'b.');
        hold on
    end
    axis equal
end
