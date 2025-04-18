function demo5()
% demo5 of im2mesh tool (convert 2d segmented image to FEM meshes)
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
    % Example 5
    % Demonstrate the case of 'tolerance = eps'
    % Polygon simplification is muted. Only knickpoints remain. Redundant 
    % vertices are removed by dpsimplify.m.
    % Note: you need to download MESH2D.
    % ---------------------------------------------------------------------
    
    % import grayscale segmented image
    im = imread('images\kumamon.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end
    
    % parameters             
    tolerance = eps;	% mute polygon simplification                  
    hmax = 500;
    mesh_kind = 'delaunay';
    grad_limit = +0.25;
   
    bounds1 = im2Bounds( im );
    bounds2 = simplifyBounds( bounds1, tolerance );

    [ node_cell, edge_cell ] = genNodeEdge( bounds2 );
    [ vert,tria,tnum ] = poly2mesh( node_cell, edge_cell, hmax, mesh_kind, grad_limit );
    
    % show result
    imshow( im,'InitialMagnification','fit' );
    plotMeshes( vert, tria, tnum );
    drawnow;
    set(figure(1),'units','normalized', ...
        'position',[.05,.35,.30,.35]);
    set(figure(2),'units','normalized', ...
        'position',[.35,.35,.30,.35]);
    % ---------------------------------------------------------------------
    
end