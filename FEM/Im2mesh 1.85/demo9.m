function demo9()
% demo9 of im2mesh tool (convert 2d segmented image to FEM meshes)
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
    % Example 9
    % Demonstrate the parameter of 'grad_limit'
    % 'grad_limit'is a scalar gradient-limit for mesh
    % Note: you need to download MESH2D.
    % ---------------------------------------------------------------------
    
    % import grayscale segmented image
    im = imread('images\b1.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end
    
    % parameters
    select_phase = [];
    tf_avoid_sharp_corner = false;  
    tolerance = 1.;
    hmax = 500;
    mesh_kind = 'delaunay';

    % case 1
    grad_limit_1 = +0.25;
    [ vert_1,tria_1,tnum_1 ] = im2mesh( im, select_phase, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit_1 );

    % case 2
    grad_limit_2 = +0.50;
    [ vert_2,tria_2,tnum_2 ] = im2mesh( im, select_phase, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit_2 );

    % show result
    plotMeshes( vert_1, tria_1, tnum_1 );
    title( ['gradient-limit of mesh = ', num2str(grad_limit_1) ])
    plotMeshes( vert_2, tria_2, tnum_2 );
    title( ['gradient-limit of mesh = ', num2str(grad_limit_2) ])
    drawnow;
    set(figure(1),'units','normalized', ...
        'position',[.05,.35,.30,.35]);
    set(figure(2),'units','normalized', ...
        'position',[.35,.35,.30,.35]);
    % ---------------------------------------------------------------------

end