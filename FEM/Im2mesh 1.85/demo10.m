function demo10()
% demo10 of im2mesh tool (convert 2d segmented image to FEM meshes)
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
    % Example 10
    % Demonstrate the parameter of 'select_phase'
    %
    % 'select_phase' is a new parameter for function im2mesh.m
    % Parameter type: vector
    % If 'select_phase' is [], all the phases will be chosen.
    % 'select_phase' is an index vector for sorted  grayscales (ascending 
    % order) in an image.
    % For example, an image with grayscales of 40, 90, 200, 240, 255. If 
    % you're interested in 40, 200, and 240, then set 'select_phase' as 
    % [1 3 4]. Those phases corresponding to grayscales of 40, 200, and 240
    % will be chosen to perform meshing.
    %
    % Note: you need to download MESH2D.
    % ---------------------------------------------------------------------
    
    % import grayscale segmented image
    im = imread('images\b1.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end
    
    % parameters
    tf_avoid_sharp_corner = false;  
    tolerance = 1.;
    hmax = 500;
    mesh_kind = 'delaunay';
    grad_limit = +0.25;

    % case 1
    select_phase_1 = [];
    [ vert_1,tria_1,tnum_1 ] = im2mesh( im, select_phase_1, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit );

    % case 2
    select_phase_2 = [3 5 7];
    [ vert_2,tria_2,tnum_2 ] = im2mesh( im, select_phase_2, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit );

    % show result
    plotMeshes( vert_1, tria_1, tnum_1 );
    plotMeshes( vert_2, tria_2, tnum_2 );
    drawnow;
    set(figure(1),'units','normalized', ...
        'position',[.05,.35,.30,.35]);
    set(figure(2),'units','normalized', ...
        'position',[.35,.35,.30,.35]);
    % ---------------------------------------------------------------------

end