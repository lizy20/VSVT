function [c,dc,g1,dg1,g2,dg2,g3,dg3]= computesenstive(F,U,opt, dwdy,dVdy,edofMat,vol_lim,cand)
        dcdx = zeros(opt.nele,1); 
        dcdwli = zeros(opt.nele,cand);
        dcdalpha = zeros(opt.nele,1);  % OBJECTIVE
        Udis= U(edofMat);
        % c=0;
        for i=1:opt.nele
             dcdx(i,1) =  -   Udis(i,:)*opt.dKe_dx( :, [8*i-7:8*i]) *Udis(i,:)';
            dcdalpha(i,1) = - Udis(i,:)*opt.dKe_dalpha(:,8*i-7:8*i)* Udis(i,:)'; 
          for j= 1:cand
            dcdwli(i,j) = - Udis(i,:)*opt.dKe_dwli(:,8*i-7:8*i,j)* Udis(i,:)';
          end  
        end

        c = F'*U;      %%========================================================== sensitive 
        dcdz = opt.P2*(dwdy.*dcdwli);    %
        dcdx= opt.P2 * dcdx; % 
        dcdalpha1 = opt.P2 * dcdalpha;  
        dc = [dcdz(:); dcdx(:) ; dcdalpha1(:) ];
        %======================================================================  %g1 real vol
       dvdx= sum(opt.V,2);  
       g1 = sum(dvdx.*opt.xPhys(:).*opt.ElemArea(:)) - vol_lim*sum(opt.ElemArea(:));
       dvdz=repmat(opt.xPhys(:), 1,cand);
       dvdz= opt.P2 *(opt.ElemArea(:) .* dVdy .* dvdz) ;
       dvdx=opt.P2 *  (opt.ElemArea(:).*dvdx);
       dvdalpha1 =zeros(opt.nele,1); 
       dg1 = [ dvdz(:); dvdx(:) ; dvdalpha1(:) ];
       %============================================================= %%g2  Z-converegen
       [g2, dg2dZ]=pnormconstraint2(opt.Z);
       dg2 = [ dg2dZ(:); 0*dvdx(:); 0*dvdalpha1(:)];

 end 


