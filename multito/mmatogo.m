-%���ĸ���top3d(30,20,10,0.4,3,1.5)

clear all;

top3d(24,24,24,0.6,3,1.5);

function top3d(nelx,nely,nelz,volfrac,penal,rmin)


% ����ѭ������
maxloop = 200;    % ����������
tolx = 0.01;      % ��ֹ�������в
displayflag = 1;  % ��ʾ�ṹ��ʾ
% ���ϵ�����
E0 = 15000;           % �������������ģ��
nu = 0.3;         % ���ɱ�

% % % % USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(0:nelx, 0:nely, nelz);                 % Coordinates
% loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
% loaddof = 3*loadnid(:) - 1;                             % DOFs

x=boundaryimport(nelx,nely,nelz,nu);

% 
loaddof=csvread('loaddof.csv');
load=csvread('Fboundary.csv');

% % USER-DEFINED SUPPORT FIXED DOFs
% [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
% fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
% fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs

% ����Ԫ��������
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);   %%x��Ӧ row�� y��Ӧcol.
F = sparse(loaddof',1,load',ndof,1);


U = zeros(ndof,1);
% [row3,~]=size(fixednid);
% for ip=1:row3
%     U(fixednid(ip)*3-2:fixednid(ip)*3,1)=fixeddis(ip*3-2:ip*3,1);
% end


%--------------
%�̶�local voi3�ǵ㣬�۲�nodeco�ֶ����롣
fixednid=[1,25,601]
fixeddof = [3*fixednid(:), 3*fixednid(:)-1, 3*fixednid(:)-2];


freedofs = setdiff(1:ndof,fixeddof);

KE = lk_H8(nu);

% ������Ԫ��ڵ���
%ik ��jkΪ�ն���
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1); %!!!!!!!!+/-1
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
% edofVec = 3*nodeids(:)+1;
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [ 3 4 5 0 1 2 ] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [ 3 4 5  0 1 2 ] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% ������



iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1; %��y��x��z
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1; %���嵥Ԫ���
                        jH(k) = e2; %��ع�������ڵ���
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));% Ȩ��ϵ��
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);


% ������ʼ��
mf= 4; %subvoxel length
m = 2; %Լ������
% x = repmat(volfrac,[nely,nelx,nelz]);
loop = 0;
change = 1;
xori=x; %volfrac Ӧ��Ϊ��ʼ����ܶ�

xPhys = x;
xmin = repmat(0.01,[nely*nelx*nelz,1]);low = xmin;
xmax = ones(nely*nelx*nelz,1);upp = xmax;
xold1 = x(:);xold2 = x(:);
n = nely*nelx*nelz;
c_MMA=1000*ones(m,1);
d=zeros(m,1);
a0=1;
a=zeros(m,1);
display_3D(xPhys);


% ��ʼ����
while change > tolx && loop < maxloop
    loop = loop+1;
    
    % ����Ԫ����
    sK = reshape(KE(:)*(xPhys(:)'.^penal*15000),24*24*nele,1);
    K = sparse(iK,jK,sK); K=(K+K')/2;    
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    
    % ����Ŀ�꺯���������ȷ���
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((xPhys.^penal*15000).*ce)));
    dc = -penal*15000*xPhys.^(penal-1).*ce;
 
    [vx,dv]=gmd(x,mf,xori);
    [gmc,dgmc]=gmdc(x,n);
          
    dv = ones(nely,nelx,nelz);
    
    % �������Ƚ��й�����ͶӰ����Ҫע����ǣ���������ѧ��֤���������Ƚ��й��˵�
    % �ȶ��ԣ�������ֲ��ʶ���������Գ��Ը��Ĺ���Ϊ�ܶȹ���
    dc(:) = H*(dc(:)./Hs);  
    dv(:) = H*(dv(:)./Hs);
    
    
    % �ƶ������߷����
    xval  = x(:);
    f0val = c;
    df0dx = dc(:);
    fval  =  [vx; gmc];   
    dfdx  = [dv(:)';dgmc(:)'];
    
    
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
    mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    % ���µ�������
    xnew     = reshape(xmma,nely,nelx,nelz);
    xPhys(:) = (H*xnew(:))./Hs;
    xold2    = xold1(:);
    xold1    = x(:);
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    % ������
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    % ����ܶ�
     if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
end
clf; display_3D(xPhys);
end

function [gm,dgm]=gmd(x,mf,xori)
    [lx,ly,lz]=size(x);
    nlx=lx/mf; nly=lx/mf;nlz=lz/mf;
    nm=nlx*nly*nlz;
    num=1;
    bb=ones(nlx,nly,nlz);
    bbori=ones(nlx,nly,nlz);
    dgm=zeros(lx,ly,lz);
    for i=1:nlx
        for j= 1:nly
            for k=1:nlz
           bb(i,j,k)= sum(sum(sum(x((i-1)*mf+1:i*mf,(j-1)*mf+1:j*mf,...
               (k-1)*mf+1:k*mf))))/(mf^3); 
           bbori(i,j,k)=sum( sum(sum(xori((i-1)*mf+1:i*mf,(j-1)*mf+1:j*mf,...
               (k-1)*mf+1:k*mf))))/(mf^3); 
           end   
        end
    end   
    kk=(norm(bb(:)-bbori(:)))^2/nm- 0.01; %% Լ����ֵ��=0.01
  
    for q=(i-1)*mf+1:i*mf
               for w= (j-1)*mf+1:j*mf
                   for e=(k-1)*mf+1:k*mf
                     dgm(q,w,e)=2*kk/(lx*ly*lz);
                   end
               end
    end
    gm=kk;
end
function [gmc,dgmc]=gmdc(x,n)
 uu=(1.-x).*(x);
 gmc= sum(uu(:))/n-0.01;
 dgmc=(1.-2*x)/n;
end



% === չʾ3DЧ��ͼ (ISO-VIEW) ===
function display_3D(rphy)

[nely,nelx,nelz] = size(rphy);
rho= rphy;    %����ԪԪ������
            
hx = 1; hy = 1; hz = 1;            % ����Ч��ͼ�ĵ�Ԫ��С
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.01)  % ����չʾ�����ܶ���ֵ
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[1-rho(j,i,k),1-rho(j,i,k),1-rho(j,i,k)]);
                hold on;
%                 [0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end