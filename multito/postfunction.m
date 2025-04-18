%========================function postdealing
function [c]= postfunction(filename,app)
 disp(' Post-process optimization ...' );
 myPrint(app,  ' Post-process optimization ...' );
load( [ pwd,'\',filename,'\final.mat' ] );
Vertices=opt.Node;
%===================================
 opt.clis = clis;
 opt.vlis = vlis;
 opt.cand =cand;
if cand >1
    1;
else
    disp(' Only one component found.' );
end

%======== Unit length for modeling,=====================================
% opt.lenx =  length of square TO unit 
% opt.lenx=opt.Lx/opt.nelx;    opt.leny=opt.lenx;
avglen=avgEdgeLength(opt.Node,opt.Element); 
opt.lenx=avglen;
opt.leny=avglen;
opt.meshlen=avglen*0.01;
%=================================================================
cp1= opt.P; cp2 =opt.P2;    
opt.P = PolyFilter(opt, max(rmin*0.5, 1.6*avglen));
opt.P2 = opt.P';
% Reduce filtering radius, to eliminate possible blank interface between materials
%Too small a value at the (rmin) position can cause checkboard problems.
%Manually adjust the value for the desired effect.
[opt, low,upp, xmin,xmax,xold1,xold2,loop] =subiterZ(conkey, opt,penal, DNN, cand, E0, Emin,gamma,loop, ci, 1, vol_lim,...
                                 low, upp, xmin, xmax, xold1, xold2, FEMparam ,0, 25,app);
opt=Zinterpo(opt,conkey);
%===========================================================================
gamma=1;
opt.P=cp1; opt.P2=cp2;
[opt, low,upp, xmin,xmax,xold1,xold2,loop] =subiterx(conkey, opt, penal, DNN, cand, E0, Emin,gamma,loop, ci, 1, vol_lim,...
                                 low,upp,xmin,xmax,xold1, xold2,  FEMparam ,0,25,app) ;
%========================================================================
opt= updateinfo(opt);
figure(1)
yyaxis left
plot(1: length(opt.clis), opt.clis(1:end), '-');
ylabel('Compliance');  
yyaxis right
plot(1: length(opt.vlis),round(opt.vlis(1:end),2),   '--');
ylabel('Volume fraction')
ylim([0.2,1]);
xlabel('iteration');
colormap(slanCM('tab20'))

if cand>1
[Color,RBM] = ConvertColors(opt.V);
figure(2); hold on;
I =reshape(Color, [],3);
patch( 'Faces',padcellrows(opt.Element), 'Vertices',Vertices,'FaceVertexCData',I ,...
     'Facecolor', 'flat', 'EdgeColor','none','HandleVisibility', 'off');
c1 = plot(nan, nan, 's',  'MarkerFaceColor', [   0.5647    0.7961    0.9843],...
      'MarkerEdgeColor', 'k', 'DisplayName', ['Aiso=', num2str(opt.ml(1))]);
c2 = plot(nan, nan, 's',  'MarkerFaceColor', [ 0.9647    0.7922    0.8980],...
      'MarkerEdgeColor', 'k', 'DisplayName', ['Aiso=', num2str(opt.ml(2))]);
c3 = plot( nan, nan, 's',  'MarkerFaceColor', [  0.5882    0.8000    0.7961],...
    'MarkerEdgeColor', 'k', 'DisplayName', ['Aiso=', num2str(opt.ml(3))]);
legend show;  axis equal;
drawnow;
end

[DH,dDHdx]=DH_fit2D(E0, Emin, opt.xPhys(:), cand,DNN);
[opt.Ke, opt.dKe_dx, opt.dKe_dalpha, opt.dKe_dwli] = elementMatVec2D(opt, gamma, DH, dDHdx, cand);
sK = reshape(opt.Ke(:),64*opt.nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2; 
U = zeros(opt.ndof,1);
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
c = F'*U;
opt.v2 = sum(opt.V,2).* opt.xPhys(:).*opt.ElemArea;
opt.v3 = sum(opt.V,2).* opt.xPhys(:);
fprintf(' It:final  Step.:%5i, H_beta.:%5.3f,  Obj.:%7.4f , Vol.:%7.3f,  finalVol.:%7.3f  \n',...
ci,1, c,  sum(opt.xPhys(:))/opt.nele,  sum(opt.v2(:))/sum(opt.ElemArea(:)) );

figure(3) %===============================================================
clist=sum(opt.V,2)<=0.35;
opt.alpha(clist) = -89;
patch(gca,'Faces',cell2mat(opt.Element),'Vertices',Vertices,'FaceVertexCData',opt.alpha,...
     'Facecolor', 'flat','EdgeColor','none');
clim([-90 90]); axis equal;  drawnow;
colormap(slanCM('Blues'))
colorbar;
figure(4) %===============================================================
opt.xPhys=opt.xPhys(:);
opt.xPhys(clist) = opt.xMin;
patch(gca,'Faces',cell2mat(opt.Element),'Vertices',Vertices,'FaceVertexCData',opt.xPhys,...
     'Facecolor', 'flat','EdgeColor','none');
 clim([0 1]); axis equal; drawnow; colorbar;
colormap( slanCM('binary'));
colorbar;

var=setdiff(who, {'app'});
save( [ pwd,'\',filename,'\final_p.mat' ] , var{:} );
end

function opt = updateinfo(opt)
opt.xPhys = opt.P*opt.x(:);  
opt.vx2 = sum(opt.V,2).*opt.xPhys(:).*opt.ElemArea;
disp(  [ 'final volume:' , num2str( sum(opt.vx2(:))/sum(opt.ElemArea))] );
end

