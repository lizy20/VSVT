function [opt, low,upp, xmin,xmax,xold1,xold2,vlis,clis,loop] = loopopt(conkey, filename, opt, DNN, cand, Param, E0, Emin,gamma, loop, ...
    ci, beta, vol_lim, low, upp, xmin, xmax, xold1, xold2, FEMparam,penal, vlis,clis,herikey,maxloop,app)

switch(conkey)
    case(1);    m=1;
    case(2);    m=2;
end
if  opt.FEMini ==0
    opt=FEMini(opt);
    opt.FEMini =1;
end
if(~exist('maxloop','var'));   maxloop=30; end

iK=FEMparam.iK ; jK= FEMparam.jK; freedofs=FEMparam.freedofs;
F =FEMparam.F; edofMat =FEMparam.edofMat ;
c_MMA=2000*ones(m,1);d=zeros(m,1);a0=1;a=zeros(m,1);
loopParam =0;  

if herikey==1
1;
else
low =  xmin;
upp = xmax;
end
n = length(xold1);
while loopParam < maxloop %&& opt.zchange>0.02
        loop = loop+1;   loopParam =loopParam+1;  
        if loopParam ==1
            opt.Y =  opt.P * opt.Z;
            opt.xPhys = opt.P* opt.x(:); % 
            opt.alpha = opt.P *opt.alpha1;
            [opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y,0.5,beta,'SIMP-H',penal);
        end
        [DH,dDHdx]=DH_fit2D(E0, Emin, opt.xPhys(:), cand,DNN);  % =====================FE-ANALYSIS AT TWO SCALES
        [opt.Ke, opt.dKe_dx, opt.dKe_dalpha, opt.dKe_dwli] = elementMatVec2D(opt, gamma, DH, dDHdx, cand);
        sK = reshape(opt.Ke(:),64*opt.nele,1);
        K = sparse(iK, jK,sK); 
        K = (K+K')/2;
        U = zeros(opt.ndof,1);
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
        c=F(freedofs,1)'*U(freedofs,1);
        opt.K=K;
        % ======================================================================
        [c,dc,g1,dg1,g2,dg2]= computesenstive(F,U,opt, dwdy,dVdy,edofMat,vol_lim,cand);
        %=======================================================================
        xval  = [opt.Z(:); opt.x(:); opt.alpha1(:)];
        f0val = c; df0dx = dc;
        switch(conkey)
        case(1)
              fval  = [g1;];  dfdx  = [ dg1'; ];
        case(2)
             fval  = [g1; g2];  dfdx  = [ dg1';dg2' ];
        end

        [xNew, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] =mmasub(m, n, loopParam, xval, xmin, xmax, xold1, xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
        opt.Z = reshape( xNew(1:cand *opt.nele), [opt.nele,cand]);
        opt.Y =  opt.P * opt.Z;   [opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y, 0.5, beta,'SIMP-H',penal);

        opt.x = xNew( 1+cand *opt.nele : (cand+1)* opt.nele);
        opt.xPhys = opt.P*opt.x(:); 

        opt.alpha1 =  xNew( 1+(1+cand)*opt.nele : (cand+2)* opt.nele) ;
        opt.alpha =  opt.P*opt.alpha1(:);
     
        % PRINT RESULTS
        opt.change= (xNew(:)-xval(:))./(xmax-xmin);      opt.zchange= max(abs( opt.change(1:cand *opt.nele,1)));
        opt.xchange= max( abs(  opt.change( 1+cand*opt.nele:(cand+1)*opt.nele,1 )));
        opt.alphachange = max(abs(  opt.change( 1+(1+cand)*opt.nele : (cand+2)* opt.nele,1 )));
        xold2  = xold1(:);   xold1  = xval(:);

        %%=============================================
        opt.vx2 = sum(opt.V,2).*opt.xPhys(:) ;
        vol_p=sum(opt.vx2(:).*opt.ElemArea)/sum(opt.ElemArea);
        fprintf(' It.:%5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , Vol.:%7.3f,  Z_ch.:%7.3f, x_ch.:%7.3f , alpha_ch.:%7.3f,   \n',...
        loop, ci, beta,c,  vol_p, opt.zchange, opt.xchange, opt.alphachange);
        myPrint(app, ...
        ' It.:%5i,  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , Vol.:%7.3f,  Z_ch.:%7.3f, x_ch.:%7.3f , alpha_ch.:%7.3f\n', ...
          loop, ci, beta, c, vol_p, opt.zchange, opt.xchange, opt.alphachange);
         
        if mod(loopParam,5)==0
            drawimage(opt,app)
        end
        if mod(loopParam,5)==0 && beta<=25 && ci== size(Param,2) -1 & loopParam>50  %10
        beta=beta+1;
        end
       
        while(app.IsPaused)
              drawnow;  % UI control pause
              pause(0.05);
       end
       if  loopParam>2 & mod(loopParam,40)==0  %40-10
          [opt, low,upp, xmin,xmax,xold1,xold2] =subiteralpha(conkey, opt, penal, DNN, cand, E0, Emin,gamma, loop, ci, ...
                 beta, vol_lim, low, upp, xmin, xmax, xold1, xold2, FEMparam,0, 20,app) ;
       end
        %%=============================================
        clis=[clis,c];
        vlis=[vlis, vol_p];
        yyaxis(app.UIAxes, 'left');
        ylabel(app.UIAxes,'Compliance')
        plot(app.UIAxes,1:length(clis), clis); 
        yyaxis(app.UIAxes, 'right');
        ylabel(app.UIAxes,'Vol_{frac}')
        plot(app.UIAxes,1:length( vlis),  vlis);
        if mod(loopParam,50)==0
            save( [ filename,'\break', num2str(loopParam),'.mat' ] );
        end
end
end

%------------------------------------------------- TABULATE SHAPE FUNCTIONS
function fem = TabShapeFnc(fem)
ElemNNode = cellfun(@length,fem.Element); % number of nodes per element

fem.ShapeFnc = cell(max(ElemNNode),1);
for nn = min(ElemNNode):max(ElemNNode)
  [W,Q] = PolyQuad(nn);
  fem.ShapeFnc{nn}.W = W;
  fem.ShapeFnc{nn}.N = zeros(nn,1,size(W,1));
  fem.ShapeFnc{nn}.dNdxi = zeros(nn,2,size(W,1));
  for q = 1:size(W,1)
    [N,dNdxi] = PolyShapeFnc(nn,Q(q,:));
    fem.ShapeFnc{nn}.N(:,:,q) = N;
    fem.ShapeFnc{nn}.dNdxi(:,:,q) = dNdxi;
  end
end
end

%------------------------------------------------ POLYGONAL SHAPE FUNCTIONs
function [N,dNdxi] = PolyShapeFnc(nn,xi)
N=zeros(nn,1); alpha=zeros(nn,1); dNdxi=zeros(nn,2); dalpha=zeros(nn,2);
sum_alpha=0.0; sum_dalpha=zeros(1,2); A=zeros(nn,1); dA=zeros(nn,2);
[p,Tri] = PolyTrnglt(nn,xi);
for i=1:nn
  sctr = Tri(i,:); pT = p(sctr,:);
  A(i) = 1/2*det([pT,ones(3,1)]);
  dA(i,1) = 1/2*(pT(3,2)-pT(2,2));
  dA(i,2) = 1/2*(pT(2,1)-pT(3,1));
end
A=[A(nn,:);A]; dA=[dA(nn,:);dA];
for i=1:nn
  alpha(i) = 1/(A(i)*A(i+1));
  dalpha(i,1) = -alpha(i)*(dA(i,1)/A(i)+dA(i+1,1)/A(i+1));
  dalpha(i,2) = -alpha(i)*(dA(i,2)/A(i)+dA(i+1,2)/A(i+1));
  sum_alpha = sum_alpha + alpha(i);
  sum_dalpha(1:2) = sum_dalpha(1:2)+dalpha(i,1:2);
end
for i=1:nn
  N(i) = alpha(i)/sum_alpha;
  dNdxi(i,1:2) = (dalpha(i,1:2)-N(i)*sum_dalpha(1:2))/sum_alpha;
end
end

%---------------------------------------------------- POLYGON TRIANGULATION
function [p,Tri] = PolyTrnglt(nn,xi)
p = [cos(2*pi*((1:nn))/nn); sin(2*pi*((1:nn))/nn)]';
p = [p; xi];
Tri = zeros(nn,3); Tri(1:nn,1)=nn+1;
Tri(1:nn,2)=1:nn; Tri(1:nn,3)=2:nn+1; Tri(nn,3)=1;
end

%----------------------------------------------------- POLYGONAL QUADRATURE
function [weight,point] = PolyQuad(nn)
[W,Q]= TriQuad;                  %integration pnts & wgts for ref. triangle
[p,Tri] = PolyTrnglt(nn,[0 0]);  %triangulate from origin
point=zeros(nn*length(W),2); weight=zeros(nn*length(W),1);
for k=1:nn
  sctr = Tri(k,:);
  for q=1:length(W)
    [N,dNds] = TriShapeFnc(Q(q,:));  %compute shape functions
    J0 = p(sctr,:)'*dNds;
    l = (k-1)*length(W) + q;
    point(l,:) = N'*p(sctr,:);
    weight(l) = det(J0)*W(q);
  end                                 
end
end

%---------------------------------------------------- TRIANGULAR QUADRATURE
function [weight,point] = TriQuad
point=[1/6,1/6;2/3,1/6;1/6,2/3]; weight=[1/6,1/6,1/6];   
end

%----------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
function [N,dNds] = TriShapeFnc(s)
N=[1-s(1)-s(2);s(1);s(2)]; dNds=[-1,-1;1,0;0,1];
end

%------------------------------------------------------------- INITIAL PLOT


function [fem] =FEMini(fem)

 fem.ElemNDof = 2*cellfun(@length,fem.Element);
 fem.i = zeros(sum(fem.ElemNDof.^2),1); 
 fem.j=fem.i; fem.k=fem.i; fem.e=fem.i;
 index = 0;
 if ~isfield(fem,'ShapeFnc'), fem=TabShapeFnc(fem); end

  for el = 1:fem.nele  
    NDof = fem.ElemNDof(el);
    eDof = reshape([2*fem.Element{el}-1; 2*fem.Element{el}],NDof,1);
    I=repmat(eDof ,1,NDof); J=I';
    fem.iK(index+1:index+NDof^2) = I(:);
    fem.jK(index+1:index+NDof^2) = J(:); 
    fem.eK(index+1:index+NDof^2) = el;
    index = index + NDof^2;
  end
    if ~isfield(fem,'ElemArea')
      fem.ElemArea = zeros(fem.nele,1);
      for el=1:fem.nele
        vx=fem.Node(fem.Element{el},1); vy=fem.Node(fem.Element{el},2);
        fem.ElemArea(el) = 0.5*sum(vx.*vy([2:end 1])-vy.*vx([2:end 1]));
      end
    end
end

