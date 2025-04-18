function plotBounds( bounds1 )
% plot bounds{i}{j}
% show its starting points and control points
    
    figure;
    hold on
    axis image off;
    
    for i = 1: length(bounds1)
        for j = 1: length(bounds1{i})

            poly = bounds1{i}{j};
            plot( poly(:,1), poly(:,2) );
            % show starting points with 'x'
            plot( poly(1,1), poly(1,2), 'kx' );
            
            % find NaN
            TFNaN = isnan( poly(:,1) );
            index = find( TFNaN ) - 1;
            % show control points with 'o'
            for k = 1: length(index)
                plot( poly(index(k),1), poly(index(k),2), 'ro' );
            end

        end
    end
    hold off
end
