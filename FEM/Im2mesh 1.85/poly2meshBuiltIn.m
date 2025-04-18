function [vert,tria,tnum] = poly2meshBuiltIn( node_cell, edge_cell, pcell, hgrad, hmax, hmin )
% poly2meshBuiltIn: generate meshes of parts defined by polygons using 
%                   matlab built-in function generateMesh
%
% input:   
% (Please check documentation of matlab built-in function generateMesh.)
%     hgrad       % Mesh growth rate
%     
%     hmax        % Target maximum mesh edge length
% 
%     hmin        % Target minimum mesh edge length
%   
% output:
%   vert(k,1:2) = [x_coordinate, y_coordinate] of k-th node 
%   tria(m,1:3) = [node_numbering_of_3_nodes] of m-th element
%   tnum(m,1) = n; means the m-th element is belong to phase n
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Jan 2022
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   
    
    [ node, edge, ~ ] = regroup( node_cell, edge_cell );
    
    % ---------------------------------------------------------------------
    % Delaunay triangulation in 2D
    DT = delaunayTriangulation(node,edge);
    % triplot(DT)

    tnodes = DT.Points';
    telements = DT.ConnectivityList';
    
    % ---------------------------------------------------------------------
    % label regionID
    regionID = zeros( size(DT.ConnectivityList,1), 1 );

    C = incenter(DT);
    for i = 1: size(C,1)
        for j = 1: length(pcell)
            if isinterior(pcell{j},C(i,1),C(i,2))   % need improve
                regionID( i ) = j;
                break;
            end
        end
    end

    regionID = regionID';
    
    % ---------------------------------------------------------------------
    % convert to geometry object in pde model
    model = createpde;
    geometryFromMesh(model,tnodes,telements,regionID);
    % pdegplot(model,'FaceLabels','on')
    
    % ---------------------------------------------------------------------
    % generate mesh using matlab built-in function
    mesh = generateMesh(model,'Hgrad',hgrad,'Hmax',hmax,'Hmin',hmin,'GeometricOrder','linear' );
    % figure
    % pdemesh(model)
    
    % ---------------------------------------------------------------------
    % obtain variable tnum
    % tnum is marker for phase
    tnum = zeros( size(mesh.Elements,2), 1 );
    count = 1;

    for i = 1: length(pcell)
        for j = 1: pcell{i}.NumRegions

            indEf = findElements(mesh,'region','face',count);
            tnum(indEf) = i;
            count = count + 1;
        end
    end
    
    % ---------------------------------------------------------------------
    % obtain variable vert, tria
    vert = mesh.Nodes';
    tria = mesh.Elements';
    % plotMeshes( vert, tria, tnum );
    % ---------------------------------------------------------------------
end

