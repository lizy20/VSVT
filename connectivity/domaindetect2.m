function [phi_2D, phi_2D_ele] = domaindetect2( points,info,le,t,ang,subidx, key,colorse )

x=points(:,1);y=points(:,2);
x0=info(1,:);y0=info(2,:);
n=length(le);
x1=x0-cosd(ang).*le/2;  x2=x0+cosd(ang).*le/2; 
y1=y0-sind(ang).*le/2; y2= y0+sind(ang).*le/2;

%x,y, column vector
%  x1,y1,x2,y2 row vector
x1=repmat(x1,[length(x),1]);
x2=repmat(x2,[length(x),1]);

y1=repmat(y1,[length(x),1]);
y2=repmat(y2,[length(x),1]);
t= repmat(t,[length(x),1]);

l = ((x2-x1).^2+(y2-y1).^2).^0.5;
sin_theta = (y2-y1)./l;
cos_theta = (x2-x1)./l;

if key==1
    
phi_r1= 1- (-1*sin_theta .*(x-x0) + cos_theta.*(y-y0) ).^2 ./ (t/2).^2;
phi_r2= 1- ((cos_theta).*(x-x0)+sin_theta.*(y-y0)).^2 ./ (l/2).^2;
phi_r= min(phi_r1,phi_r2);
[phi_2D,ind]= max(phi_r,[],2);

elseif key==2
%===========================================
trans = reshape( [cosd(ang);sind(ang)], [1,2,n]);
cood = reshape(  [x0;y0],[2,1,n]);
x_0 = reshape(pagemtimes(trans,cood),1,n);

trans = reshape( [-sind(ang);cosd(ang)], [1,2,n]);
y_0= reshape( pagemtimes(trans,cood),1,n);

trans = reshape( [cosd(ang);sind(ang)], [1,2,n]);
cood = repmat([x';y'], [1,1,n]);
xp = squeeze(pagemtimes(trans,cood));

trans = reshape( [-sind(ang);cosd(ang)],  [1,2,n]);
yp= squeeze(pagemtimes(trans,cood));

%============================================================
phi_r= 1- ((xp-x_0).^2 ./ le.^2 +(yp-y_0).^2./t.^2);
% phi_r= 1- (abs(xp-x_0)./ le + abs(yp-y_0)./t);
% phi_r= 1- max(abs(xp-x_0)./ le,abs(yp-y_0)./t);

[phi_2D,ind]= max(phi_r,[],2);
elseif key==3
1;
end
vorcolor= 2+rand(length(info),1);
% vorcolor= randi[30:10: 30+10*(length(info)-1)];
phi_2Dre=phi_2D;
phi_2D(phi_2Dre<0)=0;
phi_2D(phi_2Dre>0)=colorse(subidx(ind(phi_2Dre>0)));
phi_2D_ele = subidx(ind);
end