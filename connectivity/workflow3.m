%======================================
function workflow3(filename,app)
disp('Compute Voronoi tessellation ...');
myPrint(app,'Compute Voronoi tessellation ...');
load([pwd, '\' , filename , '\geo_data2.mat']);   
opt.kcal =1; % extra meshlen range  for block domain detection
block_height = 20;    % up -efficient // down -- memory requirement // GPU recommand =25 for 20GB
% //CPU recommand =50- 100      for 64GB-128GB
GPUkey=0;   % 1 --- GPU,0 --- CPU
% ==== LOAD & PREPROCESS ====

opt.meshlen = 0.01* avglen;     %  0.005-0.01 better qulity for numerical computation
eleid = connectfilter(points, opt);              % points including seed // no bounds

Nodex = opt.Node(:,1); Nodey = opt.Node(:,2);
Node_mx = Nodex(enidElem); Node_my = Nodey(enidElem);
opt.Xlow = (floor(min(Node_mx(:)) / avglen) ) * avglen;
opt.Xup  = (ceil( max(Node_mx(:)) / avglen) ) * avglen;
opt.Ylow = (floor(min(Node_my(:)) / avglen) ) * avglen;
opt.Yup  = (ceil( max(Node_my(:)) / avglen) ) * avglen;

[t, le, ang,  eleid, kcal, meshlen, lenx, leny, points, elelist] =   transgpu(points, eleid, opt,elelist,GPUkey,app);
[hullptsx_tol, hullptsy_tol] = meshgrid(opt.Xlow:meshlen:opt.Xup, opt.Yup:-meshlen:opt.Ylow);
result_imc_tol = zeros(size(hullptsx_tol));
result_eleid_tol = zeros(size(hullptsx_tol));
colorse = colordeal(length(points), 8);

h = waitbar(0, 'Computing Voronoi tessellation, please wait');

for j = 1:block_height:size(hullptsx_tol,1)
    if j > size(hullptsx_tol,1), break; end
    [imc, result_eleid, ~] = block_compute(opt, j,block_height, le, t, ang, kcal, meshlen, ...
                        lenx, leny, points, eleid, elelist, colorse, hullptsx_tol, hullptsy_tol);

    [row_cnt, ~] = size(imc);
    row_end = j + row_cnt - 1;
    result_imc_tol(j:row_end, :) = imc;
    result_eleid_tol(j:row_end, :) = result_eleid;
    waitbar(j/size(hullptsx_tol,1), h, sprintf('Computing Voronoi tessellation...%.2f %%', 100*j/size(hullptsx_tol,1)))
end
close(h);
figure(10);
imagesc( [opt.Xlow,opt.Xup],[opt.Yup,opt.Ylow], result_imc_tol); hold on
set (gca,'YDir',' normal ');
axis equal
patch('vertices',opt.Node,'faces',opt.enidElem,...
    'facevertexCdata',0*opt.Node(:,1),'facecolor','interp','edgecolor','k', 'facealpha',0); hold on
drawnow;
scatter(points(:,1),points(:,2),3.5,'filled', 'white'); hold on
var=setdiff(who, {'app'});
save([pwd, '\' , filename , '\geo_data3.mat'], var{:});
end


function [imc, result_eleid, hullpts] = block_compute(opt, j, block_height,le, t, ang, kcal, meshlen, ...
                                                       lenx, leny, points, eleid, elelist, colorse, ...
                                                       hullptsx_tol, hullptsy_tol)

num_rows = size(hullptsx_tol, 1);
block_end = min(j + block_height - 1, num_rows);

yfu = hullptsy_tol(j, 1);
yfl = hullptsy_tol(block_end, 1) ;
yl_ele = max(yfl, opt.Ylow);
yu_ele = min(yfu, opt.Yup);
xl_ele = opt.Xlow;
xu_ele = opt.Xup;

mask = points(:,1) >= xl_ele- kcal*opt.lenx & points(:,1) <= xu_ele + kcal*opt.lenx  & ...
       points(:,2) >= yl_ele - kcal*opt.leny & points(:,2) <= yu_ele +  kcal*opt.leny;

subpoints = points(mask, :);
subt      = t(mask);
suble     = le(mask);
subang    = ang(mask);
subidx    = find(mask);

[xnode, ynode] = meshgrid(xl_ele:meshlen:xu_ele, yu_ele:-meshlen:yl_ele);
info = [xnode(:), ynode(:)];
[phi_2D, phi_2D_ele] = domaindetect2( info,subpoints', suble, subt, subang, subidx, 2, colorse);
if ~isempty(phi_2D)
imc = reshape(phi_2D, size(xnode));
result_eleid = reshape(phi_2D_ele, size(xnode));
else
imc =[];
result_eleid=[];
end
hullpts = [xnode(:), ynode(:)];
end

function eleid = connectfilter(refcood,opt)
Mod= opt;
enidElem= opt.enidElem;
i=1;
eleid=ones(size(refcood,1),1);
for i = 1:length(enidElem)
lis = enidElem(i,:);
in = inpolygon(refcood(:,1),refcood(:,2),opt.Node(lis,1),opt.Node(lis,2));
eleid(in)=i;
end
end


function [t, le, ang,  eleid, kcal, meshlen, lenx, leny, points, elelist] ...
    = transgpu(points,eleid,opt,elelist,gpukey,app)
%% gpu array interface for 
np= size(points,1);
Pdfilter=DistPntSets(points,points, 0.4*opt.lenx );
t=10* ones( 1,np );
le=reshape( opt.aniso(eleid),1, length(t));
le=le.*t;
ang=reshape( opt.alpha(eleid) , 1, length(t));
le=  (Pdfilter* le')'; 
ang=  (Pdfilter* ang')';
kcal =opt.kcal; meshlen=opt.meshlen;
lenx=opt.lenx ; leny=opt.leny ;
dig =40;
if gpukey ==1
try
    g = gpuDevice();  % 
    if ~isempty(g)
        disp(['Using GPU: ', g.Name]);
        myPrint(app, ['Using GPU: ', g.Name]);
    ang= gpuArray(ang);  le=gpuArray(le); t=gpuArray(t);
    eleid = gpuArray(eleid); 
    kcal=gpuArray(opt.kcal);
    meshlen = gpuArray(opt.meshlen);
    lenx=gpuArray(opt.lenx) ; leny=gpuArray(opt.leny) ;
    elelist=gpuArray(elelist);
    else
         disp('No usable GPU detected.Using default CPU');
         myPrint(app, 'No usable GPU detected.Using default CPU');
    end
catch
    1;
end
end
end

function  colorlist = colordeal(x,bin)
colorf= linspace(2,10,bin); 
colorlist = kron( ones(floor(x/bin),1), colorf');
if mod(x,bin)~=0
colorlist= [colorlist; colorf(1: mod(x,bin))' ];
end
end

function [P] = DistPntSets(PS1,PS2,R)
SElem=[];
d = cell(size(PS1,1),1);
for el = 1:size(PS1,1)       %Compute the distance information
    dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    [I,J] = find(dist<=R);   %Find the indices for distances less that R
    [~,K] = setdiff(I,SElem);
    d{el} = [I(K),J(K)+(el-1),dist(I(K))];
end
d = cell2mat(d);   
P = sparse(d(:,1),d(:,2),1-d(:,3)/R);    %Assemble the filtering matrix
P = spdiags(1./sum(P,2),0,size(PS1,1),size(PS1,1))*P;
end

function interpolated_surface_with_mask(x, y, z, varargin)
p = inputParser;
addOptional(p, 'InterpMethod', 'natural');  % 'linear', 'cubic' 
addOptional(p, 'GridRes', 100);         
addOptional(p, 'ViewMode', '3D');           
parse(p, varargin{:});
interpMethod = p.Results.InterpMethod;
gridRes = p.Results.GridRes;
viewMode = p.Results.ViewMode;

xq = linspace(min(x), max(x), gridRes);
yq = linspace(min(y), max(y), gridRes);
[Xq, Yq] = meshgrid(xq, yq);
x = gather(x); y = gather(y); z = gather(z);
Xq = gather(Xq);Yq = gather(Yq);
Zq = griddata(x, y, z, Xq, Yq, interpMethod);
shp = alphaShape(x, y, inf); 
inMask = inShape(shp, Xq, Yq);

Zq(~inMask) = NaN;

surf(Xq, Yq, Zq, 'EdgeColor', 'none');
colorbar;
axis equal;
switch lower(viewMode)
    case '2d'
        view(2);
    case '3d'
        view(3);
end
end



