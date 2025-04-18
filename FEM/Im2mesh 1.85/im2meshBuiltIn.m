function [ vert, tria, tnum ] = im2meshBuiltIn( im, tf_avoid_sharp_corner, tolerance, hgrad, hmax, hmin )
% im2meshBuiltIn: convert grayscale segmented image to triangular mesh  
%                 using matlab built-in function generateMesh
%
% input
%     im        % grayscale segmented image
%
%     tf_avoid_sharp_corner   % For function getCtrlPnts
%                             % Whether to avoid sharp corner when 
%                             % simplifying polygon.
%                             % Value: true or false
%                             % If true, two extra control points
%                             % will be added around one original 
%                             % control point to avoid sharp corner 
%                             % when simplifying polygon.
%                             % Sharp corner in some cases will make 
%                             % poly2mesh not able to converge.
%                         
%     tolerance   % For funtion simplifyBounds
%                 % Tolerance for polygon simplification.
%                 % Check Douglas-Peucker algorithm.
%                 % If u don't need to simplify, try tolerance = eps.
%                 % If the value of tolerance is too large, some 
%                 % polygons will become line segment after 
%                 % simplification, and these polygons will be 
%                 % deleted by function delZeroAreaPoly.
%     
% (Please check documentation of matlab built-in function generateMesh for 
%  parameter hgrad, hmax, and hmin. 
%  https://www.mathworks.com/help/pde/ug/pde.pdemodel.generatemesh.html)
%
%     hgrad       % For funtion poly2meshBuiltIn
%                 % Mesh growth rate
%     
%     hmax        % For funtion poly2meshBuiltIn
%                 % Target maximum mesh edge length
% 
%     hmin        % For funtion poly2meshBuiltIn
%                 % Target minimum mesh edge length
%   
% output:
%   vert(k,1:2) = [x_coordinate, y_coordinate] of k-th node 
%   tria(m,1:3) = [node_numbering_of_3_nodes] of m-th element
%   tnum(m,1) = n; means the m-th element is belong to phase n
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Jan 2020
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   
   
    bounds1 = im2Bounds( im );
    bounds2 = getCtrlPnts( bounds1, tf_avoid_sharp_corner );
    % plotBounds( bounds2 );
    
    bounds3 = simplifyBounds( bounds2, tolerance );
    bounds3 = delZeroAreaPoly( bounds3 );
    % plotBounds( bounds3 );

    % clear up redundant vertices
    % only control points and knick-points will remain
    bounds4 = getCtrlPnts( bounds3, false );
    bounds4 = simplifyBounds( bounds4, eps );
    
    [ node_cell, edge_cell ] = genNodeEdge( bounds4 );
    pcell = bound2polyshape(bounds4);
    % generate mesh
    [vert,tria,tnum] = poly2meshBuiltIn( node_cell, edge_cell, pcell, hgrad, hmax, hmin );
    % plotMeshes( vert, tria, tnum );
end

