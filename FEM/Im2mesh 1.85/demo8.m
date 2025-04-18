function demo8()
% demo8 of im2mesh tool (convert 2d segmented image to FEM meshes)
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
    % Example 8
    % Demonstrate the parameter of 'tf_avoid_sharp_corner'
    %
    % Whether to avoid sharp corner when simplifying polygon.
    % Value: true or false
    % If true, two extra control points will be added around one original 
    % control point to avoid sharp corner when simplifying polygon.
    % Inspiration: sharp corner in some cases will make function poly2mesh 
    % not able to converge.
    % Note: you need to download MESH2D.
    % ---------------------------------------------------------------------
    
    % import grayscale segmented image
    im = imread('images\c1.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end

    % parameters
    select_phase = [];
    tolerance = 1.1;
    hmax = 500;
    mesh_kind = 'delaunay';
    grad_limit = +0.25;
    
    % case 1
    tf_avoid_sharp_corner = 0;
    % use evalc() to suppress command window output of im2mesh
    evalc([ '[ vert_1,tria_1,tnum_1 ] = im2mesh( im, select_phase, ' ...
            'tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit )' 
            ]);

    % case 2
    tf_avoid_sharp_corner = 1;
    evalc([ '[ vert_2,tria_2,tnum_2 ] = im2mesh( im, select_phase, ' ...
            'tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit )' 
            ]);

    % show result
    plotMeshes( vert_1, tria_1, tnum_1 );
    title( ['tf\_avoid\_sharp\_corner = 0, |TRIA|=', ...
        num2str(size(tria_1,1)) ])
    plotMeshes( vert_2, tria_2, tnum_2 );
    title( ['tf\_avoid\_sharp\_corner = 1, |TRIA|=', ...
        num2str(size(tria_2,1)) ])
    drawnow;
    set(figure(1),'units','normalized', ...
        'position',[.05,.35,.30,.35]);
    set(figure(2),'units','normalized', ...
        'position',[.35,.35,.30,.35]);
    % ---------------------------------------------------------------------
    
end