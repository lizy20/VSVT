function [] = workflow4(filename,app)
% Voronoi numerical pattern --> truss branches info
myPrint(app,'Compute truss geometry ...');
load([pwd, '\' , filename , '\geo_data3.mat']);
% vertex= zeros(length(points),100);
%==============================================normal
M=length(points);
% h=waitbar(0,'Computing Voronoi truss info, please wait...');
hullpts_tol= [hullptsx_tol(:), hullptsy_tol(:)];
idxhull =[1:length(hullpts_tol)]';
re_filter = result_eleid_tol(:);
barinfo=[];
%%==========================================================
% % % %%% quick boudary smooth 
% % % % % load('MBBbonudspoly.mat','c5','c6'); c5,c6 boundary nodes
% % % % % nbounds_poly= subtract(polyshape(c5.Position), polyshape(c6.Position));
% % % % % nbounds_poly.simplify();
% % % % % plot(nbounds_poly)
% % % % % [inB,onB] =  isinterior( nbounds_poly , hullptsx_tol(:), hullptsy_tol(:) );
%----------------------------------------------------------------------------------------------------------
pg = polyshape(nbounds(:,1), nbounds(:,2));
inB = isinterior(pg, hullptsx_tol(:), hullptsy_tol(:));
re_filter(~(inB),:)= 0;
%============================================================
result_eleid_tol = reshape(re_filter, size(result_eleid_tol,1),size(result_eleid_tol,2) ) ;
barinfo = runTrussExtractionParallel(hullpts_tol, M,result_eleid_tol, opt, idxhull,4);
%%%%%%======================================================
% close(h);
% figure(12);
% line([barinfo(1:200,1), barinfo(1:200,3)]',[barinfo(1:200,2), barinfo(1:200,4)]' );
var=setdiff(who, {'app'});
save([pwd, '\' , filename , '\geo_data4.mat'],  var{:});
end


function [ K1,subbarinfo]= calhull(points,idf)
K= convhull(gather(points));
convex_hull_points = points(K, :);
boundary_lines = [convex_hull_points(1:end-1, :), convex_hull_points(2:end, :)];
% plot( points(K,1), points(K,2), 'b-'); hold on;
threshold =5; 
slopes = atan2d( boundary_lines(:,4)- boundary_lines(:,2),   boundary_lines(:,3)- boundary_lines(:,1) );
% merged_boundary_lines = [];
Sdiff = abs(slopes(2:end)-slopes(1:end-1));
idx =[1:length(Sdiff)]'+1;
K1=K;  K1(idx(Sdiff<threshold),:)=[]; 
% merged_boundary_lines =[points(K1(1:end-1),:),  points(K1(2:end),:)   ];
subbarinfo=[ points(K1(1:end-1),1), points(K1(1:end-1),2), points(K1(2:end),1), points(K1(2:end),2) ];
% plot( points(K1,1), points(K1,2), 'k-'); hold on; %%  draw sketch
K1= [idf(K1)', zeros(1, 100-length(K1))];
end

function barinfo = quick_bar_extract(points, alpha, simplifyTol)
% quick_bar_extract
% Extract approximate straight-line bars along the boundary of a point set.
%
% Inputs:
%   points       - [N x 2] coordinates of 2D scattered points
%   alpha        - boundary sensitivity parameter (0 < alpha <= 1)
%                  alpha ~1: convex-like; alpha ~0: tight boundary
%   simplifyTol  - simplification tolerance for reducing vertices (optional)
%
% Output:
%   barinfo      - [M x 4] straight-line bars as [x1 y1 x2 y2]
    if nargin < 2
        alpha = 1;  % default to convex-like
    end
    if nargin < 3
        simplifyTol = 0.5;  % default tolerance
    end
    % Step 1: Compute approximate boundary
    K = boundary(points, alpha);
    B = points(K, :);
    % Step 2: Convert to polyshape and simplify
    pg = polyshape(reducepoly([B(:,1), B(:,2)], simplifyTol));  % reduce noise
    % Step 3: Extract simplified vertices
    V = pg.Vertices;
    if isempty(V)
        barinfo = []; return;
    end
    % Step 4: Convert to line segments
    barinfo = [V(1:end-1,:), V(2:end,:)];  % [x1 y1 x2 y2]
    % Optional: close loop
    if norm(V(end,:) - V(1,:)) > 1e-6
        barinfo(end+1,:) = [V(end,:), V(1,:)];  % close last edge
    end
end
                               
function barinfo = runTrussExtractionParallel(hullpts_tol, M,result_eleid_tol, opt, idxhull, numWorkers)
    if nargin < 6
        numWorkers = 4;  
    end
    barinfoCell = cell(M, 1);
    
    delete(gcp('nocreate'));
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= numWorkers
        delete(gcp('nocreate'));
        parpool(numWorkers);
    end
   WaitMessage = parfor_wait_ui(M,'ReportInterval', ceil(M/20),'Waitbar', true, 'Title', 'Computing truss info ... '); 
    parfor i = 1:M
        c5 = result_eleid_tol(:) == i;
        temp = [];
        if any(c5)
            try
                % temp = quick_bar_extract(hullpts_tol(c5,:),0.95, opt.meshlen*1); 
                [~, temp ] =calhull( hullpts_tol(c5,:) ,idxhull(c5) );  
            catch
                warning('Failed at region %d', i);
            end
        end
        barinfoCell{i} = temp;

      
        WaitMessage.Send;   
        drawnow;


    end
    barinfo = vertcat(barinfoCell{:});
    WaitMessage.Destroy;
end
