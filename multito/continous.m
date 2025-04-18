function []= continous(filename, Macro_struct,vol_lim,cand ,ml ,Fkey,rmin,app)
% element ： convex ， tri or quai
% elemAt： node must in clockwise or anti-clockwise
warning OFF

%==================Basic refered material propoerty=========================
%================Continuation on material interpolation parameters==========
E0 = 1;  Emin = 1e-9; nu = 0.3;
Param = [1,2,3,3,3,3;   
            0,0.25,0.5,0.75,1,1]; %========= =
conkey=1;
xi=0.5; 
mkdir([pwd, '\', filename]);
%=================Mesh data infor :Typical Rectangle element domain or nonstructure meshinp======
opt.Lx = Macro_struct(1); opt.Ly = Macro_struct(2);
opt.nelx   = Macro_struct(3); opt.nely  = Macro_struct(4);
opt.iniVol    = Macro_struct(5); 
opt.Elex = opt.Lx/opt.nelx; opt.Eley = opt.Ly/opt.nely;
opt.nele = opt.nelx*opt.nely;  
opt.ndof = 2*(opt.nelx+1)*(opt.nely+1);
opt.FEMini =0;
opt.nonreg=0; %default 0 --> regular mesh
% ===================================================================
if ~ismember(Fkey, [5,6,7]) 
[x,y]=meshgrid(  0:opt.Elex:opt.Lx,opt.Ly:-opt.Eley:0);
opt.Node=[x(:),y(:)];
nodenrs = reshape(1:(opt.nely+1)*(opt.nelx+1),1+opt.nely,1+opt.nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,opt.nele,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*opt.nely+[2 3 0 1] -2 -1],opt.nele,1);
opt.Element = num2cell(edofMat(:,2:2:end)/2 , 2);
end

switch(Fkey)
% %=========cantilever========================
 case(1)
[load_x, load_y] = meshgrid( opt.nelx, (opt.nely)/2);
loadnid  = round(load_x*(opt.nely+1)+(opt.nely+1-load_y));
F = sparse( 2*loadnid(:), 1, -1, 2*(opt.nelx+1)*(opt.nely+1),1);
U = zeros(opt.ndof,1);
[fixed_x, fixed_y] = meshgrid(0, 0:opt.nely);
fixednid  = fixed_x*(opt.nely+1)+(opt.nely+1-fixed_y);
fixeddofs = [2*fixednid(:); 2*fixednid(:)-1];
freedofs  = setdiff(1:opt.ndof,fixeddofs);
% % %==========================================
%=============MBB============================
 case(2)
[load_x, load_y] = meshgrid( 0, opt.nely);
loadnid  = load_x*(opt.nely+1)+(opt.nely+1-load_y);
F = sparse( 2*loadnid(:), 1, -1,  2*(opt.nelx+1)*(opt.nely+1),1);
U = zeros(opt.ndof,1);
[fixed_x, fixed_y] = meshgrid(opt.nelx, 0);
fixednid1  = fixed_x*(opt.nely+1)+(opt.nely+1-fixed_y);
[fixed_x, fixed_y] = meshgrid(0,0:opt.nely);
fixednid2  = fixed_x*(opt.nely+1)+(opt.nely+1-fixed_y);
fixednid=[fixednid1;fixednid2];
fixeddofs = [2*fixednid1(:);fixednid2(:)*2-1];
freedofs  = setdiff(1:opt.ndof,fixeddofs);
%======================================================
%=============uni-press==================================
   case{3,4}
[load_x, load_y] = meshgrid(0:opt.nelx, opt.nely);
loadnid  = load_x*(opt.nely+1)+(opt.nely+1-load_y);
F = sparse( 2*loadnid(:), 1, -1/(opt.nelx+1),  2*(opt.nelx+1)*(opt.nely+1),1);
U = zeros(opt.ndof,1);
[fixed_x, fixed_y] = meshgrid(0:opt.nelx,0);
fixednid  = fixed_x*(opt.nely+1)+(opt.nely+1-fixed_y);
fixeddofs = [2*fixednid(:); fixednid(:)*2-1];
freedofs  = setdiff(1:opt.ndof,fixeddofs);

%======================================================
%==================User-defined mesh info===================
case(5)
opt.nonreg=1;
Meshinfo=load([pwd,'\meshdata\FEMmeshdata.mat']);
opt.Node = Meshinfo.node;   
opt.nele = size(Meshinfo.element,1);
opt.ndof =size(Meshinfo.node,1)*2;  

edofMat = kron(Meshinfo.element, [2,2]);      edofMat(:, 1:2:end)= edofMat(:, 1:2:end)-1;
opt.Element = num2cell(Meshinfo.element,2);   
fixednid = Meshinfo.fixednid;           fixeddofs =Meshinfo.fixeddofs;
loadnid = Meshinfo.loadnid;             loaddofs = Meshinfo.loadndofs;

freedofs  = setdiff(1:opt.ndof,fixeddofs);
F = sparse(   loaddofs, 1, Meshinfo.Fvalue,  opt.ndof,1);
U = zeros(opt.ndof,1);
%===============Demo-1========================
case(6)
opt.nonreg=1;
Meshinfo=load([pwd,'\meshdata\FEMmeshdata_1.mat']);
opt.Node = Meshinfo.node;   
opt.nele = size(Meshinfo.element,1);
opt.ndof =size(Meshinfo.node,1)*2;  

edofMat = kron(Meshinfo.element, [2,2]);      edofMat(:, 1:2:end)= edofMat(:, 1:2:end)-1;
opt.Element = num2cell(Meshinfo.element,2);   
fixednid = Meshinfo.fixednid;           fixeddofs =Meshinfo.fixeddofs;
loadnid = Meshinfo.loadnid;             loaddofs = Meshinfo.loadndofs;

freedofs  = setdiff(1:opt.ndof,fixeddofs);
F = sparse(   loaddofs, 1, Meshinfo.Fvalue,  opt.ndof,1);
U = zeros(opt.ndof,1);
%==================Demo-2=======================
case(7)
opt.nonreg=1;
Meshinfo=load([pwd,'\meshdata\FEMmeshdata_2.mat']);
opt.Node = Meshinfo.node;   
opt.nele = size(Meshinfo.element,1);
opt.ndof =size(Meshinfo.node,1)*2;  

edofMat = kron(Meshinfo.element, [2,2]);      edofMat(:, 1:2:end)= edofMat(:, 1:2:end)-1;
opt.Element = num2cell(Meshinfo.element,2);   
fixednid = Meshinfo.fixednid;           fixeddofs =Meshinfo.fixeddofs;
loadnid = Meshinfo.loadnid;             loaddofs = Meshinfo.loadndofs;

freedofs  = setdiff(1:opt.ndof,fixeddofs);
F = sparse(   loaddofs, 1, Meshinfo.Fvalue,  opt.ndof,1);
U = zeros(opt.ndof,1);
end

figure(10)
title('Mesh overview')
patch( 'vertices',opt.Node,'faces', padcellrows(opt.Element),...
    'facevertexCdata',0*opt.Node(:,1),'facecolor','interp','edgecolor','k', 'facealpha',0);
hold on 
axis equal
scatter(opt.Node(fixednid,1), opt.Node(fixednid,2), 50, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
scatter(opt.Node( [loadnid],1),opt.Node( [loadnid],2), 50, 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 1.5);
legend('Mesh Element','Boundary Nodes', 'Loading Nodes');
drawnow;

iK = reshape(kron(edofMat,ones(8,1))', 64*opt.nele, 1 );
jK = reshape(kron(edofMat,ones(1,8))', 64*opt.nele, 1 );
opt.cand =cand;
opt.ml=ml;
Z=0.8*ones(opt.nele,cand);
xPhys=1*ones(opt.nele,1);


% INITIALIZE ITERATION
opt.alpha1 = zeros(opt.nele,1);
opt.x = xPhys;
opt.Z=Z;
opt.change = 1;
opt.zMin=0; opt.zMax=1;
opt.xMin= Macro_struct(6); 
opt.xMax=Macro_struct(7);
opt.alphaMin= -89.5; opt.alphaMax=90;
if Fkey==4
    opt.alphaMin= 89; opt.alphaMax=90; %% Extra alpha constraint corresponding case in DOI:  doi.org/10.1016/j.cma.2024.117378
end
MaxIter= [20,20,20,20,150,0];   %   [50,50,50,50,100,0];  
vlis=[];clis=[];
loop = 0; 

switch(conkey)
    case(1)
        m=1;
    case(2)
        m=2;
    case(3)
        m=3;
end

xmin = [ opt.zMin * ones(opt.nele*cand,1); 
              opt.xMin * ones(opt.nele,1);
              opt.alphaMin * ones(opt.nele,1);];
xmax = [opt.zMax * ones(opt.nele*cand,1); 
              opt.xMax * ones(opt.nele,1);
              opt.alphaMax *ones(opt.nele,1);];
low = xmin;upp = xmax;n = length(xmax);
%===========================================================================
for ki= 1:cand
    for did= 1:4
        net=load([pwd,'\multito\grad\NNplanestressangD_0.2_0.9_', num2str(ml(ki)),   '_',num2str(did)  ,'.mat']);
        DNN{ki,did}=net; 
    end
end
FEMparam.iK =iK  ;  FEMparam.jK = jK; 
FEMparam.freedofs=freedofs;
FEMparam.F=F;  FEMparam.edofMat =edofMat; 
xold1 = [opt.Z(:); opt.x(:); opt.alpha1(:)];
xold2 = xold1;
%=================================================================
beta=1; gamma =0.1;
for ci =1:size(Param,2)%          5 % 
      penal = Param(1,ci);
      gamma = Param(2,ci); %--------------------------------------------------------------------gamma set
      maxloop = MaxIter(ci);
      opt.P = PolyFilter(opt,rmin);
      opt.P2 = opt.P';
      disp( ['current p: ', num2str(penal),...
        ', current gamma: ', num2str(gamma)]);
    opt.zchange = 1;
    opt.zMin=0; opt.zMax=1;  
    opt.ZPRMove=0.05; opt.ZPREta=0.2;
    [opt, low,upp, xmin,xmax,xold1,xold2,vlis,clis,loop] = loopopt(conkey,filename, opt, DNN, cand,Param, E0, Emin,gamma, loop, ...
        ci, beta, vol_lim, low, upp, xmin, xmax, xold1, xold2, FEMparam, penal, vlis,clis, 1,maxloop,app);
end
var=setdiff(who, {'app'});
save( [ pwd,'\',filename,'\final.mat' ] , var{:} );






