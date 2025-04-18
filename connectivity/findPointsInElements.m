function points_in_elements = findPointsInElements(midpoints, enidElem, Node)
numElements = size(enidElem, 1);
points_in_elements = cell(numElements, 1);
for i = 1:numElements
    idx = enidElem(i,:);
    verts = Node(idx, :);
    in = inpolygon(midpoints(:,1), midpoints(:,2), verts(:,1), verts(:,2));
    points_in_elements{i} = find(in);
end
end
