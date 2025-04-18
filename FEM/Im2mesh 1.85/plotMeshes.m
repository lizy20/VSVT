function plotMeshes( vert, tria, tnum )
% show meshes

    figure;
    hold on; 
    axis image off;

    tvalue = unique( tnum );
    num_phase = length( tvalue );
    
    % setup color
    if num_phase == 1
        col = 0.98;
    elseif num_phase > 1
        col = 0.3: 0.68/(num_phase-1): 0.98;
    else
        error("num_phase < 1")
    end
    
    for i = 1: num_phase
        phasecode = tvalue(i);
        patch('faces',tria( tnum==phasecode, 1:3 ),'vertices',vert, ...
        'facecolor',[ col(i), col(i), col(i) ], ...
        'edgecolor',[.1,.1,.1]);
    end
    hold off
    
%     drawnow;
%     set(figure(1),'units','normalized', ...
%         'position',[.05,.50,.30,.35]) ;
end

