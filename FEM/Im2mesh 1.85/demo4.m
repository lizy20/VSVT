function demo4()
% demo4 of im2mesh tool (convert 2d segmented image to FEM meshes)
% 
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Oct 2020
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   

    % download MESH2D from www.mathworks.com/matlabcentral/fileexchange/25
    % 555-mesh2d-delaunay-based-unstructured-mesh-generation, and then add
    % the folder (mesh2d-master) to your path using addpath function
    addpath('mesh2d-master')
    
    % ---------------------------------------------------------------------
    % Example 4
    % Demonstrate what is inside function im2mesh
    % Note: you need to download MESH2D.
    % ---------------------------------------------------------------------
    
    % import grayscale segmented image
    im = imread('images\kumamon.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end
    
    % parameters
    tf_avoid_sharp_corner = false;  % Whether to avoid sharp corner
    tolerance = 1.;                 % Tolerance for polygon simplification
    hmax = 500;                     % Maximum mesh-size
    mesh_kind = 'delaunay';         % Meshing algorithm
    grad_limit = +0.25;             % Scalar gradient-limit for mesh
    
    
    % image to polygon boundary
    bounds1 = im2Bounds( im );
    bounds2 = getCtrlPnts( bounds1, tf_avoid_sharp_corner );
    
    % simplify polygon boundary
    bounds3 = simplifyBounds( bounds2, tolerance );
    bounds3 = delZeroAreaPoly( bounds3 );

    % clear up redundant vertices
    % only control points and knick-points will remain
    bounds4 = getCtrlPnts( bounds3, false );
    bounds4 = simplifyBounds( bounds4, eps );
    
    % generate mesh
    [ node_cell, edge_cell ] = genNodeEdge( bounds4 );
    [ vert,tria,tnum ] = poly2mesh( node_cell, edge_cell, hmax, mesh_kind, grad_limit );
    % plotMeshes( vert, tria, tnum );
    
    % show result
    imshow( im,'InitialMagnification','fit' );
    plotBounds( bounds2 );
    plotBounds( bounds4 );
    plotMeshes( vert, tria, tnum );
    drawnow;
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]);
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]);
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]);
    set(figure(4),'units','normalized', ...
        'position',[.35,.05,.30,.35]);
    % ---------------------------------------------------------------------

end