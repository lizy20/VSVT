function p = bound2polyshape(bounds)
    p = cell( length(bounds), 1 );  % p is a cell vector for polyshape
    
    % boundary cell to polyshape cell
    for i = 1: length(bounds)
        x_temp = [];
        y_temp = [];
        for j = 1: length(bounds{i})
            x_temp = [ x_temp; NaN; bounds{i}{j}(:,1) ];
            y_temp = [ y_temp; NaN; bounds{i}{j}(:,2) ];
        end
        p{i} = polyshape( x_temp, y_temp );
        p{i} = simplify( p{i} );
    end
end

