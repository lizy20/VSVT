function demo1()
% demo1 of im2mesh tool (convert 2d segmented image to FEM meshes)
% 
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Jan 2022
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   

    % ---------------------------------------------------------------------
    % Example 1
    % Demonstrate function im2meshBuiltIn, which use matlab built-in 
    % function generateMesh() to generate mesh from geometry
    % Note: you don't need to download MESH2D for this example.
    % ---------------------------------------------------------------------
    % Note: 
    % Im2mesh program identifies different phases in a image by their 
    % grayscales. Different grayscales correspond to different phases. If 
    % you have 4 level of grayscales in a image, the resulted meshes will 
    % contain 4 phases.
    % You need to do pre-processing to your image before you put the image 
    % into Im2mesh program. For example, convert to 8-bit grayscale image, 
    % noise removal, and image segmentation (e.g., Otsu's method). What 
    % Im2mesh needs is a segmented grayscale image. 
    % ---------------------------------------------------------------------
    
    % import grayscale segmented image
    im = imread('images\kumamon.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end
    
    % parameters
    tf_avoid_sharp_corner = false;  % Whether to avoid sharp corner
    tolerance = 1.;                 % Tolerance for polygon simplification
    hgrad = 1.5;                    % Mesh growth rate
    hmax = 100;                     % Target maximum mesh edge length
    hmin = 0.01;                    % Target minimum mesh edge length
    % the detailed explanation of these parameters is in im2meshBuiltIn.m
    
    % convert grayscale segmented image to triangular mesh using matlab 
    % built-in function generateMesh
    [ vert,tria,tnum ] = im2meshBuiltIn( im, tf_avoid_sharp_corner, tolerance, hgrad, hmax, hmin );
    % plotMeshes( vert, tria, tnum );
    
    imshow( im,'InitialMagnification','fit' );
    plotMeshes( vert, tria, tnum );
    drawnow;
    set(figure(1),'units','normalized', ...
        'position',[.05,.35,.30,.35]);
    set(figure(2),'units','normalized', ...
        'position',[.35,.35,.30,.35]);
    % ---------------------------------------------------------------------
    
end