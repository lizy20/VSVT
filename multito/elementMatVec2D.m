%% SUB FUNCTION: elementMatVec2D

 function [Ke,dKe_dx,dKe_dalpha,dKe_dwli]= elementMatVec2D(opt,gamma,DH_t,dDH_tdx,m)

% a=opt.Elex/2; b= opt.Eley/2; 
Wb = opt.W; 
nele= opt.nele;

Ke = zeros(8,8*nele);     
dKe_dx=zeros(8,8*nele);
dKe_dalpha=zeros(8,8*nele);  
dKe_dwli=zeros(8,8*nele,m);
% multimaterial interpolation function
Wc = 1-gamma*Wb;

for cc= 1:m
      Wj(:,cc) = prod( Wc(:, setdiff( [1:m],cc) ) ,2);
      Wo(:,cc) = Wb(:,cc).*Wj(:,cc);
      for cc2 =1:m
          Wr(:,cc,cc2) = prod( Wc(:, setdiff(1:m,[cc,cc2])),2);
      end
end

alpha2= -opt.alpha(:);
 
for lo= 1:nele 
    Dl=zeros(3,3); 
    dDl_dwli=zeros(3,3,m); 
    dDl_dx=zeros(3,3);
    dDl_dalpha=zeros(3,3);
    %=================================================
 % alpha stand slope, anticlockwise. 
% here DH' is the aimed ,  known;
% cal DH  under stand global coodr DH =M* D' *MT 

    R = [ cosd( alpha2(lo)) , sind( alpha2(lo)) ,0;
             -sind( alpha2(lo)) , cosd( alpha2(lo)) ,0;
             0,0,1 ];

    M=[ R(1,1)^2, R(1,2)^2,  2*R(1,1)*R(1,2);
          R(1,2)^2,   R(2,2)^2, 2*R(2,2)*R(2,1); 
         R(1,1)*R(2,1),  R(1,2)*R(2,2), cosd(2* alpha2(lo)) ] ;

    dMdalpha= -1* [ -sind(2* alpha2(lo)) , sin( 2*alpha2(lo)) , 2*cosd(2*alpha2(lo) );
        sin( 2*alpha2(lo)) , -sind(2* alpha2(lo)) , -2*cosd(2*alpha2(lo) );
        - cosd(2*alpha2(lo) ), cosd(2*alpha2(lo) ),  -2* sind(2*alpha2(lo)) ];

    %%===================Dl========================== 
    for cc=1:m 
        DHi=  reshape(DH_t(:,lo,cc),3,3);
        Dl = Dl+ Wo(lo,cc)*M* DHi* M';  
    end    
    %%===================dDldxl==========================
    for cc=1:m
        dDHidx=  reshape(dDH_tdx(:,lo,cc),3,3);
        dDl_dx = dDl_dx +  Wo(lo,cc)*M* dDHidx*M';  
    end
    %%===================dDldwli===========================
    for i = 1:m   % dDl_dwli is (3,3,m)
        DHi =  reshape(DH_t(:,lo,i),3,3);
        dDl_dwli(:,:,i)= Wj(lo,i)*M*DHi*M' ;
        for p = setdiff( [1:m], i)
             DHp= reshape(DH_t(:,lo,p),3,3);
             dDl_dwli(:,:,i) = dDl_dwli(:,:,i) - gamma * Wb(lo,p) * Wr(lo,i,p) *M* DHp*M'; %-===========================- ,+ 更改
        end
    end
    %%===================dDlidalpha==========================
    for cc=1:m
        DHi=  reshape(DH_t(:,lo,cc),3,3);
        dDli_dalpha= dMdalpha * DHi * M' + M * DHi *(dMdalpha)'; 
        dDl_dalpha = dDl_dalpha+ Wo(lo,cc)*dDli_dalpha;        
    end

nn =  length( opt.Element{lo} );
w = opt.ShapeFnc{nn}.W;
enode = opt.Element{lo};
for q=1:length(w)
    %==================================================================
      dNdxi = opt.ShapeFnc{nn}.dNdxi(:,:,q);
      J = opt.Node(enode,:)'*dNdxi; 
      dNdx = dNdxi/J;
      B = zeros(3,2*nn);
      B(1,1:2:2*nn) = dNdx(:,1)'; 
      B(2,2:2:2*nn) = dNdx(:,2)'; 
      B(3,1:2:2*nn) = dNdx(:,2)'; 
      B(3,2:2:2*nn) = dNdx(:,1)';

    % stiffness matrix
    Ke(:,[lo*8-7:lo*8]) = Ke(:, [lo*8-7:lo*8]) + B'*Dl*B*w(q)*det(J);  
    dKe_dalpha(:,[lo*8-7 :lo*8]) = dKe_dalpha(:, [lo*8-7:lo*8])  + w(q) *det(J)*B'* dDli_dalpha *B;
    dKe_dx(:, [lo*8-7:lo*8]) =  dKe_dx(:, [lo*8-7:lo*8])+ w(q)*det(J)*B'* dDl_dx *B;
            for p=1:m
                dKe_dwli(:, [lo*8-7:lo*8], p)=  dKe_dwli(:, [lo*8-7:lo*8],p)+ w(q)*det(J)*B'* dDl_dwli(:,:,p) *B;           
            end
end

end
end


