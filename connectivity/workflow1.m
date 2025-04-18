
function []= workflow1(filename,app)
 % samet: defualt same truss thickness for TO unit dehomogenazation
 % samet: defualt same seed density for TO unit dehomogenazation
 %  Nc  // default_t 

 myPrint(app,'Compute dehomogenazation parameters...');
 load( [ pwd,'\',filename,'\final_p.mat' ] );

%===============Passive element  /Material domain setting===============
opt.ElemAct1 = 1:opt.nele;
opt.ElemAct1(opt.v3<0.3) = [];
opt.nele= length(opt.ElemAct1);
opt.V=opt.V(opt.ElemAct1,:); opt.Z=opt.Z(opt.ElemAct1,:);
opt.xPhys = opt.xPhys(opt.ElemAct1);
opt.alpha = opt.alpha(opt.ElemAct1);

%================================================
opt.N=zeros(opt.nele,1); opt.t=zeros(opt.nele,1);  
opt.alphachange=0;
anisolist=opt.ml;
for i= 1:length(anisolist)
    if length(opt.ml)==1
    readdata = load( ['NNplanestressang3V_' , num2str(anisolist(1)), '.mat'] );
    RV{i} = readdata.net1;
    else
    readdata = load( ['NNplanestressang3V_' , num2str(anisolist(i)), '.mat'] );
    RV{i} = readdata.net1;
    end
end
% [opt,elelist] = samet(RV,anisolist,opt);
[opt,elelist] = sameN(RV,anisolist,opt);
recheck = ones(1,length(elelist));

for ci= 1:length(anisolist)
i =anisolist(ci);  
n1= opt.aniso(elelist); n2 =opt.N(elelist)'; n3= opt.t(elelist)'; n4 =opt.xPhys(:); n4= n4(elelist)';
if any(n1==i)
recheck(1, n1==i) = RV{ci}([n2(n1==i); n3(n1==i)]) - n4(n1==i) ;
end
end
[~,idc]= find(abs(recheck) >0.05*vol_lim);

if size(idc,2)==0;
    disp('success dehomo-parameter computation ~~');
else
    disp('warning dehomo-parameter error existing !');
    recheck(1,idc)
    n4(1,idc)
    n2(1,idc)
    n3(1,idc)
end
var=setdiff(who, {'app'});
save([pwd, '\' , filename , '\geo_data1.mat'], var{:});
end


% %================================================ same t
function [opt,elelist]= samet(RV,anisolist,opt)
Nmin=35;  Nmax=90;  
ck = sum ( repmat(anisolist, opt.nele,1).* opt.V,2);
passive = ck<0.5;    [~ ,Mi]= max(opt.V,[],2);  
opt.aniso=anisolist(Mi) ;  opt.aniso(passive)=0;
opt.anisonetidx = Mi; opt.anisonetidx(passive) = [];
elelist  = [1:opt.nele];           elelist( passive)=[];
x0= (Nmin+Nmax)/2;
xPhys=opt.xPhys(:); 
options = optimoptions('fminunc','Algorithm','quasi-newton');
M= length(elelist);
h=waitbar(0,'please wait');
% nl=patterncal(opt);
   nl=ones(M,1);
for idx=1:M
    str=['computing Voronoi seeds...',num2str(idx/M*100),'%'];    waitbar(idx/M,h,str)
    i= elelist(idx);
    pp = opt.anisonetidx(idx);
     if nl(i)==0
       default_t=0.022;   
    else
       default_t=0.03;  
    end
    count=0;
    Nc=round( fsolve(@(x)calvol(x, default_t, xPhys(i), RV{pp} ), x0,options) );
    while ( Nc<Nmin || Nc>Nmax ) && count<40
            if Nc < Nmin
                    default_t=default_t-0.001;          %0.005
            elseif Nc > Nmax
                    default_t=default_t+0.001;
            end
              Nc=round( fsolve(@(x) calvol(x, default_t, xPhys(i), RV{pp} ), x0,options) );
            count=count+1; 
    end
    if count >30
    warning('seed over threshold');
    end
    opt.N(i)= Nc; opt.t(i)= default_t;
end    
close(h);
end

%==========================================     same N
function [opt,elelist]=sameN(RV,anisolist,opt) 
% for better compliance/volume computation accuracy， recommand Nc>45;
% for better efficency, use Nc=20~30
tmin=0.01; tmax=0.12;
ck = sum ( repmat(anisolist, opt.nele,1).* opt.V,2);
passive = ck<0.5;    [~ ,Mi]= max(opt.V,[],2);  
opt.aniso=anisolist(Mi) ;    opt.aniso( passive)=0;
opt.anisonetidx = Mi; opt.anisonetidx(passive) = [];
elelist  = [1:opt.nele];        elelist( passive)=[];
x0= (tmin+tmax)/2;
xPhys=opt.xPhys(:); options = optimoptions('fminunc','Algorithm','quasi-newton');
M= length(elelist);
h=waitbar(0,'please wait');
for idx=1:M
 str=['computing ..',num2str(idx/M*100),'%'];    waitbar(idx/M,h,str)
    i= elelist(idx);
    pp = opt.anisonetidx(idx);
    Nc= 35;  
    count=0;     
   t= fsolve(@(t)calvol(Nc, t, xPhys(i), RV{pp} ), x0,options) ;
  while ( t<tmin || t>tmax ) && Nc>=25 && Nc<=105
            if  t <  tmin
                    Nc=Nc-1;
            elseif t > tmax
                    Nc=Nc+1;
            end
            t= fsolve(@(t)calvol(Nc, t, xPhys(i), RV{pp} ), x0,options) ;
            count=count+1; 
  end
    if count >90
    warning('seed over threshold,please adjusting threshold, or the volume error may occur');
    end
    opt.N(i)= Nc; opt.t(i)= t;
end
close(h);
end
%====================================================================
function  vol= calvol(x,y,c0,RV )
vol=abs ( RV([x,y]')-c0' ) ;
end


