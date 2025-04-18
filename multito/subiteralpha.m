function [opt, low,upp, xmin,xmax,xold1,xold2] =subiteralpha(conkey, opt, penal, DNN, cand, E0, Emin,gamma, loop, ...
    ci, beta, vol_lim, low, upp, xmin, xmax, xold1, xold2, FEMparam, herikey,maxloop,app ) 


m=1;
iK=FEMparam.iK ; jK= FEMparam.jK; freedofs=FEMparam.freedofs;
F =FEMparam.F; edofMat = FEMparam.edofMat;   nele =opt.nele;

c_MMA=2000*ones(m,1);  d=zeros(m,1);   a0=1;  a=zeros(m,1);  loopParam =0;  
if(~exist('maxloop','var'))
   maxloop=25;  % default 25 subloops
end

xold1c= xold1(nele*(cand+1)+1: nele*(cand+2));
xold2c =xold2( nele*(cand+1)+1: nele*(cand+2));

xminc =xmin(nele*(cand+1)+1: nele*(cand+2));
xmaxc =xmax(nele*(cand+1)+1: nele*(cand+2));

if herikey==1
lowc=low(nele*(cand+1)+1: nele*(cand+2)); % efficient search- old gradient direction
uppc=upp(nele*(cand+1)+1: nele*(cand+2));
else
lowc =  xminc; % new search-complex problem, new gradient direction
uppc = xmaxc;
end

n= length(xold1c);
[opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y,0.5,beta,'SIMP-H',penal);
while loopParam < maxloop 
% FE-ANALYSIS AT TWO SCALES
loopParam =loopParam+1;  
[DH,dDHdx]=DH_fit2D(E0, Emin, opt.xPhys(:), cand,DNN);
[opt.Ke, opt.dKe_dx, opt.dKe_dalpha, opt.dKe_dwli] = elementMatVec2D(opt, gamma, DH, dDHdx, cand);
sK = reshape(opt.Ke(:),64*opt.nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U= zeros(opt.ndof,1);
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);

% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS   
[c,dc,g1,dg1,g2,dg2] = computesenstive(F,U,opt, dwdy,dVdy,edofMat,vol_lim,cand);
%%=============== sensitive ==================== 
f0val = c;   df0dx = dc(nele*(cand+1)+1: nele*(cand+2));
fval  = -1;  dfdx  = [ 0*df0dx'  ]; % --------indepedent;
xval=[opt.alpha1(:)];
[xNew, ~, ~, ~, ~, ~, ~, ~, ~, lowc,uppc] =mmasub(m, n, loopParam, xval, xminc, xmaxc, xold1c, xold2c,f0val,df0dx, fval, dfdx, lowc, uppc, a0, a, c_MMA,d);
opt.alpha1 =  xNew ;
opt.alpha = opt.P*opt.alpha1(:);
% PRINT RESULTS
opt.alphachange = max(abs(opt.alpha1(:)-xold1c)) ;
xold2c    = xold1c(:);      xold1c    = xval(:);
fprintf(' It.:%5i,   alphasubloop: %5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , alpha_ch.:%7.3f \n',...
loop, loopParam, ci, beta,c, opt.alphachange);
 myPrint(app, ...
' It.:%5i,   alphasubloop: %5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , alpha_ch.:%7.3f \n',...
loop, loopParam, ci, beta,c, opt.alphachange);
if mod(loopParam,5)==0
       drawimage(opt,app)
end
while(app.IsPaused)
      drawnow;  % UI control pause
      pause(0.05);
end
end
xold1(nele*(cand+1)+1: nele*(cand+2)) = xold1c;
xold2( nele*(cand+1)+1: nele*(cand+2)) = xold2c;
xmin(nele*(cand+1)+1: nele*(cand+2)) = xminc ;
xmax(nele*(cand+1)+1: nele*(cand+2)) = xmaxc;
low(nele*(cand+1)+1: nele*(cand+2)) = lowc;
upp(nele*(cand+1)+1: nele*(cand+2)) = uppc;

end