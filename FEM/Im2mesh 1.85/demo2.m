function demo2()
% demo2 of im2mesh tool (convert 2d segmented image to FEM meshes)
% 
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Jan 2022
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
    % Example 2
    % Demonstrate function im2mesh, which use MESH2D to generate mesh from 
    % geometry. 
    % Note: you need to download MESH2D.
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
    select_phase = [];  % select certain phases in image for meshing
                        % Parameter type: vector
                        % If 'select_phase' is [], all the phases will be
                        % chosen.
                        % 'select_phase' is an index vector for sorted 
                        % grayscales (ascending order) in an image.
                        % For example, an image with grayscales of 40, 90,
                        % 200, 240, 255. If u're interested in 40, 200, and
                        % 240, then set 'select_phase' as [1 3 4]. Those 
                        % phases corresponding to grayscales of 40, 200, 
                        % and 240 will be chosen to perform meshing.
                        
    tf_avoid_sharp_corner = false;  % For function getCtrlPnts
                                    % Whether to avoid sharp corner when 
                                    % simplifying polygon.
                                    % Value: true or false
                                    % If true, two extra control points
                                    % will be added around one original 
                                    % control point to avoid sharp corner 
                                    % when simplifying polygon.
                                    % Sharp corner in some cases will make 
                                    % poly2mesh not able to converge.
                        
    tolerance = 1.0;    % Key parameter affecting meshing result
                        % For funtion simplifyBounds
                        % Tolerance for polygon simplification.
                        % Check Douglas-Peucker algorithm (https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm)
                        % If u don't need to simplify, try tolerance = eps.
                        % If the value of tolerance is too large, some 
                        % polygons will become line segment after 
                        % simplification, and these polygons will be 
                        % deleted by function delZeroAreaPoly.
    
    hmax = 500;             % For funtion poly2mesh
                            % Maximum mesh-size
    
    mesh_kind = 'delaunay'; % For funtion poly2mesh
                            % Meshing algorithm
                            % Value: 'delaunay' or 'delfront' 
                            % "standard" Delaunay-refinement or "Frontal-Delaunay" technique
    
    grad_limit = +0.25;     % For funtion poly2mesh
                            % Scalar gradient-limit for mesh                    

    
    % convert grayscale segmented image to triangular mesh
    [ vert,tria,tnum ] = im2mesh( im, select_phase, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit );
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