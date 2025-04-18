%%60-120 *pi/4 = 45-90
function [] = workflow2(filename,app)
myPrint(app,'Compute Voronoi seeding ...');

load([pwd, '\' , filename , '\geo_data1.mat']);
disp('Begin seeding...')
if length(opt.ml)==1
    opt.aniso= opt.ml(1)* ones(size(opt.aniso,1) ,size(opt.aniso,2));
end
opt.meshlen=opt.lenx*0.01;
enidElem=cell2mat(opt.Element(opt.ElemAct1,:));
pd=opt.Node(unique(enidElem(:)),:);

bounds_K =boundary(pd,0.91);   
% 0.91 - adjustable        % need check to ensure compact, accurate boundary.
nbounds=pd(bounds_K,:);
% plot(nbounds(:,1)' ,nbounds(:,2)' ); hold on
% scatter(pd(:,1),pd(:,2)); hold on
% axis equal
nbounds2= [min(nbounds(:,1)), max(nbounds(:,1));
    min(nbounds(:,2)), max(nbounds(:,2)) ];

%====================================
opt.enidElem=enidElem;
min_r = 0.04*opt.lenx; %--> increase, for higher efficiency but lower seeding accuracy.
fh=@(xc) var_dense(xc, opt);
points = poissonDisc2( fh, 'bounds', nbounds2, 'min_r', min_r ,'verbose',1);
while length(points)<=10
points = poissonDisc2( fh, 'bounds', nbounds2, 'min_r', min_r ,'verbose',1); 
end
%double repeat codes for correct ; single the function handle may not work well;--matlab verison error
%======
% [in,on] = inpolygon(points(:,1),points(:,2),nbounds(:,1),nbounds(:,2));
% in = in|on;
% points =points(in,:);
pg = polyshape(nbounds(:,1), nbounds(:,2));
in = isinterior(pg, points(:,1),points(:,2));
points =points(in,:);
%=================================
figure(8);
title('Voronoi Seeding Distribution')
plot(nbounds(:,1)' ,nbounds(:,2)' ); hold on;
scatter(points(:,1),points(:,2),1.5,'filled');grid on;axis equal;
axis equal
%===================================
var=setdiff(who, {'app'});
save( [pwd,'\' , filename , '\geo_data2.mat'], var{:} );
end

% =====================structure build shape change===============================
function disr= var_dense(refcood,opt)
Mod= opt;
N=Mod.N(:);
enidElem= opt.enidElem;
i=1;
eleid=ones(size(refcood,1),1);
for i = 1:length(enidElem)
lis = enidElem(i,:);
in = inpolygon(refcood(:,1),refcood(:,2),opt.Node(lis,1),opt.Node(lis,2));
eleid(in)=i;
end
S=Mod.ElemArea(opt.ElemAct1);
disr = opt.lenx/sqrt(2)./ (N(eleid)).^0.5; 
end

