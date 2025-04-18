function demo3()
% demo3 of im2mesh tool (convert 2d segmented image to FEM meshes)
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
    % Example 3
    % Demonstrate how to export mesh as inp, bdf, and .node/.ele file
    % Note: you need to download MESH2D.
    % ---------------------------------------------------------------------
    % step 1: obtain mesh using function im2mesh
    
    % import grayscale segmented image
    im = imread('images\kumamon.tif');
    if size(im,3) == 3;  im = rgb2gray( im ); end
    
    % parameters
    select_phase = [];              % select phases for meshing
    tf_avoid_sharp_corner = false;  % Whether to avoid sharp corner
    tolerance = 1.;                 % Tolerance for polygon simplification
    hmax = 500;                     % Maximum mesh-size
    mesh_kind = 'delaunay';         % Meshing algorithm
    grad_limit = +0.25;             % Scalar gradient-limit for mesh
    
    % convert grayscale segmented image to triangular mesh
    [ vert,tria,tnum ] = im2mesh( im, select_phase, tf_avoid_sharp_corner, tolerance, hmax, mesh_kind, grad_limit );
    % plotMeshes( vert, tria, tnum );

    % ---------------------------------------------------------------------
    % step 2: write inp, bdf, and .node/.ele file 
    
    % parameters for mesh export
    dx = 1; dy = 1;     % scale of your imgage
                        % dx - column direction, dy - row direction
                        % e.g. scale of your imgage is 0.11 mm/pixel, try
                        %      dx = 0.11; and dy = 0.11;
                        
    ele_order = 1;      % for getNodeEle, this variable is either 1 or 2
                        % 1 - linear / first-order element
                        % 2 - quadratic / second-order element
                        % note: printBdf only support linear element
                        
    ele_type = 'CPS3';  % element type, for printInp
    
    precision_nodecoor = 8; % precision of node coordinates, for printInp
                            % e.g. precision_nodecoor=4, dx=1 -> 0.5000 
                            %      precision_nodecoor=3, dx=0.111, -> 0.055
    
    % scale node coordinates
    vert( :, 1 ) = vert( :, 1 ) * dx;
    vert( :, 2 ) = vert( :, 2 ) * dy;
    
    % get node coordinares and elements from mesh
    [ nodecoor_list, nodecoor_cell, ele_cell ] = getNodeEle( vert, tria, ...
                                                        tnum, ele_order );
    
    % write file
    % inp file (Abaqus)
    % print as multiple parts
    printInp_multiPart( nodecoor_cell, ele_cell, ele_type, precision_nodecoor );
    % print as multiple sections
    printInp_multiSect( nodecoor_list, ele_cell, ele_type, precision_nodecoor );
    
    % bdf file (Nastran bulk data)
    printBdf( nodecoor_list, ele_cell, precision_nodecoor );
    
    % .node/.ele file 
    % haven't been tested
    printTria( vert, tria, tnum, precision_nodecoor );  
    
    % ---------------------------------------------------------------------

end