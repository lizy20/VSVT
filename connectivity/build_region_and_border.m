function [mainShape, innerBorder] = build_region_and_border(Node, Element, thickness, simplifyTol)
% build_region_and_border
% Constructs a polygonal domain from mesh elements, optionally simplifies the boundary,
% and generates an inner border by buffering outward and intersecting inward.
%
% Inputs:
%   Node         - [n x 2] node coordinates
%   Element      - [m x k] element connectivity (triangle or quad)
%   thickness    - border thickness (positive scalar; 0 means no border)
%   simplifyTol  - simplification tolerance (0 disables; default = 0.01)
%
% Outputs:
%   mainShape    - polyshape of the region (may include holes)
%   innerBorder  - polyshape of internal border area (empty if thickness = 0)
%
% Author: Zeyang Li

if nargin < 4
    simplifyTol = 0.01;
end
% Validate thickness
if ~isnumeric(thickness) || ~isscalar(thickness) || ~isfinite(thickness)
    warning('Invalid thickness. Setting to 0.');
    thickness = 0;
end

% Step 1: Create individual polyshape patches from elements
nElem = size(Element, 1);
pgs(nElem) = polyshape();  % preallocate
for i = 1:nElem
    id = Element(i,:);
    coords = Node(id, :);
    if any(isnan(coords(:)))
        pgs(i) = polyshape();  % invalid element
    else
        pgs(i) = polyshape(coords(:,1), coords(:,2));
    end
end

% Step 2: Union all patches into one full region
polyIn = union(pgs);
nRegions = polyIn.NumRegions+polyIn.NumHoles;

% Step 3: Simplify zigzag boundary (if requested)
if simplifyTol > 0
allParts =polyshape();
for k = 1:nRegions
try
    [xrings, yrings] = boundary(polyIn,k);
    part  = polyshape(reducepoly([xrings, yrings], simplifyTol));  % reduce each ring
    allParts(k) = part;
catch
    warning('Ring %d in region %d simplification failed.', i, k);
end
end
% Step 2: Merge all simplified parts
mainShape =subtract(allParts(1), union([allParts(2:end)]));
else
mainShape  = polyIn;
end
% plot(mainShape);
% Step 4: Return if no thickness
if thickness <= 0
    innerBorder = [polyshape()];  % placeholder
    return;
else
% Step 5: Buffer outward and intersect to get inner border
try
 outBuffer =polyshape();
for k = 1:nRegions
    [xx,yy]= boundary(mainShape,k);
    outBuffer(k) = polybuffer( [xx,yy],  'lines', 2 * thickness);
end   
     % plot(outBuffer);
    innerBorder = union(outBuffer);
    innerBorder = simplify(innerBorder);  % clean up
    % plot(innerBorder); hold on
    innerBorder = intersect( innerBorder, mainShape);
    %plot(innerBorder); 
catch
    warning('Buffer or intersection failed. Returning empty innerBorder.');
    innerBorder = [polyshape()];
end
end

end

