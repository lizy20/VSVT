function [opt, low,upp, xmin,xmax,xold1,xold2,loop] =subiterZ(conkey, opt, penal, DNN, cand, E0, Emin,gamma,loop, ci, beta, vol_lim,...
     low,upp,xmin,xmax,xold1, xold2,FEMparam,herikey,maxloop,app) 

switch(conkey)
    case(1);    m=1;
    case(2);    m=2;
end

freedofs=FEMparam.freedofs; edofMat =FEMparam.edofMat;
nele =opt.nele;

c_MMA=2000*ones(m,1);
d=zeros(m,1);
a0=1;
a=zeros(m,1);
loopParam = 0;  
if(~exist('maxloop','var'))
   maxloop=30;  % 如果未出现该变量，则对其进行赋值
end

xold1c = xold1(1:nele*cand);  xold2c = xold2(1:nele*cand);
xminc = xmin(1:nele*cand); xmaxc =  xmax(1:nele*cand);
if herikey==1
lowc =  low(1:nele*cand);
uppc = upp(1:nele*cand);
else
lowc =  xminc;
uppc = xmaxc;
end
n = length(xold1c);

[DH,dDHdx]=DH_fit2D(E0, Emin, opt.xPhys(:), cand,DNN);
while loopParam < maxloop 
% FE-ANALYSIS AT TWO SCALES
loopParam =loopParam+1;  loop=loop+1;
[opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y,0.5,beta,'SIMP-H',penal);
[opt.Ke, opt.dKe_dx, opt.dKe_dalpha, opt.dKe_dwli] = elementMatVec2D(opt, gamma, DH, dDHdx, cand);
sK = reshape(opt.Ke(:),64*opt.nele,1);
K = sparse(FEMparam.iK,FEMparam.jK,sK); K = (K+K')/2;
U= zeros(opt.ndof,1);
U(freedofs,:) = K(freedofs,freedofs) \ FEMparam.F(freedofs,:);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS 
 [c,dc,g1,dg1,g2,dg2]= computesenstive(FEMparam.F,U,opt, dwdy,dVdy,edofMat,vol_lim,cand);

xval  = [opt.Z(:); ];
f0val = c;  df0dx = dc(1:nele*cand);
dg1dx= dg1(1:nele*cand);   
dg2dx= dg2(1:nele*cand);

switch(conkey)
case(1)
      fval  = [g1;];  dfdx  = [ dg1dx'; ];
case(2)
     fval  = [g1; g2];  dfdx  = [  dg1dx';dg2dx' ];
end

[xNew, ~, ~, ~, ~, ~, ~, ~, ~, lowc,uppc] =mmasub( m, n, loopParam, xval, xminc, xmaxc, xold1c, xold2c,f0val,df0dx,fval,dfdx,lowc,uppc,a0,a,c_MMA,d);
opt.Z = reshape( xNew, [opt.nele,cand]);
opt.Y =  opt.P * opt.Z;
[opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y,0.5,beta,'SIMP-H',penal);

% PRINT RESULTS
opt.change= max(abs( (xNew(:)-xval(:))./(xmaxc-xminc)) ) ;
xold2c    =  xold1c(:);
xold1c    =  xval(:);
opt.vx2 = sum(opt.V,2).*opt.xPhys(:);
fprintf(' It.:%5i, Zsubloop: %5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , Vol.:%7.3f,  finalVol:%7.3f,  Z_ch.:%7.3f  \n',...
loop, loopParam, ci, beta,c, sum(opt.xPhys(:))/opt.nele,   sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea),  opt.change);
 myPrint(app, ...
' It.:%5i, Zsubloop: %5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , Vol.:%7.3f,  finalVol:%7.3f,  Z_ch.:%7.3f  \n',...
loop, loopParam, ci, beta,c, sum(opt.xPhys(:))/opt.nele,   sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea),  opt.change);
 
if mod(loopParam,1)==0
      drawimage(opt,app)
end
opt.clis= [opt.clis, c];
opt.vlis= [opt.vlis, sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea)];
yyaxis(app.UIAxes, 'left');
plot(app.UIAxes,1:length(opt.clis), opt.clis); 
yyaxis(app.UIAxes, 'right');
plot(app.UIAxes,1:length(opt.vlis),  opt.vlis);

end
xold1(1:nele*cand)=xold1c;
xold2(1:nele*cand) =xold2c ;
xmin(1:nele*cand) =xminc ;
xmax(1:nele*cand) = xmaxc;
low(1:nele*cand)  =  lowc;
upp(1:nele*cand) =  uppc;
end