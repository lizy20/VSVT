function avgLen = avgEdgeLength(Nodes, Elements)
    Elements=cell2mat(Elements);
    E = sort([Elements(:,[1,end]), Elements(:,1:end-1)]')';  % All edges (closed polygons)
    Edges = unique(sort([E(:,1), E(:,2)], 2), 'rows');       % emov rep;
    n = size(Edges,1);
    lens = zeros(n,1);
    for i = 1:n
        p1 = Nodes(Edges(i,1), :);
        p2 = Nodes(Edges(i,2), :);
        lens(i) = norm(p1 - p2);
    end
    avgLen = mean(lens);
end