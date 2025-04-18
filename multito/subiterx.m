function [opt, low,upp, xmin,xmax,xold1,xold2,loop] =subiterx(conkey, opt, penal, DNN, cand, E0, Emin,gamma,loop, ci, beta, vol_lim,...
     low,upp,xmin,xmax,xold1, xold2,FEMparam,herikey,maxloop,app)   

iK=FEMparam.iK ; jK= FEMparam.jK; freedofs=FEMparam.freedofs;
F =FEMparam.F; edofMat =FEMparam.edofMat ;
switch(conkey)
    case(1);    m=1;
    case(2);    m=1;
end
c_MMA=2000*ones(m,1);  d=zeros(m,1);    a0=1;   a=zeros(m,1);   loopParam =0;  
if(~exist('maxloop','var'))
   maxloop=30;  
end
nele =opt.nele;
xold1c = opt.x(:);  xold2c = opt.x(:);
xminc = xmin(1+ nele*cand : (cand+1)*nele);xmaxc =  xmax(1+ nele*cand : (cand+1)*nele);

if herikey==1
lowc =  low(1+ nele*cand : (cand+1)*nele);
uppc = upp(1+ nele*cand : (cand+1)*nele);
else
lowc =  xminc;
uppc = xmaxc;
end
n = length(xold1c);
[opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y,0.5,beta,'SIMP-H',penal);

 while loopParam < maxloop %&& opt.zchange>0.02
        loop = loop+1;   loopParam =loopParam+1;  
        [DH,dDHdx]=DH_fit2D(E0, Emin, opt.xPhys(:), cand,DNN);  % ============FE-ANALYSIS AT TWO SCALES
        [opt.Ke, opt.dKe_dx, opt.dKe_dalpha, opt.dKe_dwli] = elementMatVec2D(opt, gamma, DH, dDHdx, cand);
        sK = reshape(opt.Ke(:),64*opt.nele,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U = zeros(opt.ndof,1);
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
        %==================================================
        [c,dc,g1,dg1,g2,dg2]= computesenstive(F,U,opt, dwdy,dVdy,edofMat,vol_lim,cand);
        %==================================================
        xval  = [opt.x(:) ];
        dc =dc(1+ nele*cand : (cand+1)*nele);
        dg1 =dg1(1+ nele*cand : (cand+1)*nele);
        f0val = c; df0dx = dc;
      
        switch(conkey)
        case(1)
              fval  = g1;  dfdx  =  dg1'; 
        case(2)
             fval  = [g1; g2];  dfdx  = [ dg1';dg2' ];
        end

        [xNew, ~, ~, ~, ~, ~, ~, ~, ~, lowc,uppc] =mmasub(m, n, loopParam, xval, xminc, xmaxc, xold1c, xold2c,f0val,df0dx,fval,dfdx,lowc,uppc,a0,a,c_MMA,d);
        opt.x = xNew;
        opt.xPhys = opt.P*opt.x(:);  

        % PRINT RESULTS
        opt.xchange= max( xNew- xold1c);

        xold2c = xold1c(:);   xold1c= xval(:);
        
        opt.vx2 = sum(opt.V,2).*opt.xPhys(:);
        fprintf(' It.:%5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f ,  finalVol.:%7.3f,  x_ch.:%7.3f    \n',...
        loop, ci, beta,c,  sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea), opt.xchange);
        myPrint(app, ...
            ' It.:%5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f ,  finalVol.:%7.3f,  x_ch.:%7.3f    \n',...
              loop, ci, beta,c,  sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea), opt.xchange);

        opt.clis=[opt.clis,c];
        opt.vlis=[opt.vlis,sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea)];
        opt.vlis= [opt.vlis, sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea)];
        yyaxis(app.UIAxes, 'left');
        plot(app.UIAxes,1:length(opt.clis), opt.clis); 
        yyaxis(app.UIAxes, 'right');
        plot(app.UIAxes,1:length(opt.vlis),  opt.vlis);
        
end

xold1(1+ nele*cand : (cand+1)*nele) = xold1c;
xold2(1+ nele*cand : (cand+1)*nele) = xold2c ;
xmin(1+ nele*cand : (cand+1)*nele) = xminc ;
xmax(1+ nele*cand : (cand+1)*nele) = xmaxc;
low(1+ nele*cand : (cand+1)*nele)  =  lowc;
upp(1+ nele*cand : (cand+1)*nele) =  uppc;
end
