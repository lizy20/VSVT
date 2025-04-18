function  workflow5(filename,app)
% plot porous structure geometry

load([pwd, '\' , filename , '\geo_data4.mat']);
myPrint(app,'Plot porous structure geometry ... ');

% Step 1:
nBars = size(barinfo,1);
barthickness = zeros(nBars,1);
mid_barthickness = [(barinfo(:,1)+barinfo(:,3))/2, (barinfo(:,2)+barinfo(:,4))/2];
pl = findPointsInElements(mid_barthickness, enidElem, opt.Node);
for ci = 1:length(pl)
    idn = pl{ci};
    barthickness(idn,:) = opt.t(ci) * avglen;
end

% % Step 2: truss analysis
xlist = zeros(nBars, 5);
ylist = zeros(nBars, 5);
x1 = double(gather(barinfo(:,1)));
y1 = double(gather(barinfo(:,2)));
x2 = double(gather(barinfo(:,3)));
y2 = double(gather(barinfo(:,4)));
thickness  = double(gather(barthickness(:))); %
dx = x2 - x1; dy = y2 - y1;
len = hypot(dx, dy);
valid = len >= 1e-8 & thickness> 0 & len < avglen*5;   % get rid of abnormal 

x1 = x1(valid); y1 = y1(valid);
x2 = x2(valid); y2 = y2(valid);
thickness = thickness(valid);
dx = dx(valid); dy = dy(valid); len = len(valid);
n = numel(x1);
nx = -dy ./ len;
ny =  dx ./ len;
r  = thickness/ 2;
v1x = x1 + r .* nx; v1y = y1 + r .* ny;
v2x = x1 - r .* nx; v2y = y1 - r .* ny;
v3x = x2 - r .* nx; v3y = y2 - r .* ny;
v4x = x2 + r .* nx; v4y = y2 + r .* ny;
xlist = [v1x, v2x, v3x, v4x, v1x];
ylist = [v1y, v2y, v3y, v4y, v1y];

% figure(1); hold on;axis equal
% for i=1:100
% plot(xlist(i,:), ylist(i,:)); 
% end
barShapes = build_bar_shapes_parfor2([x1,y1,x2,y2,thickness]);
barShapes=simplify(barShapes);
% Step 3: draw poly
% barShapes = build_bar_shapes_parfor(xlist, ylist);
% Step 4:  union poly
fprintf('Union poly...\n');
tic
polyout =  batch_union_all(barShapes,4);
toc
var=setdiff(who, {'app'});
save([pwd, '\' , filename , '\geo_data5.mat'], var{:});

end



