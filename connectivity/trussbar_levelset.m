function [phi_2D] = trussbar_levelset(points,info,t)
x = points(:,1);y=points(:,2);
x1 = info(1,:); y1=info(2,:);
x2 = info(3,:);y2=info(4,:);
t = repmat(t,[1,length(x1)]);
%x,y, column vector
%  x1,y1,x2,y2 row vector
x1=repmat(x1,[length(x),1]);
x2=repmat(x2,[length(x),1]);

y1=repmat(y1,[length(x),1]);
y2=repmat(y2,[length(x),1]);
t= repmat(t,[length(x),1]);

l= ((x2-x1).^2+(y2-y1).^2).^0.5;
sin_theta = (y2-y1)./l;
cos_theta = (x2-x1)./l;
x0=(x1+x2)/2; y0= (y1+y2)/2;
phi_c1=-((x-x1).^2+(y-y1).^2);
phi_c2=  -((x-x2).^2+(y-y2).^2);
phi_r1= 1-(-1*sin_theta .*(x-x0) + cos_theta.*(y-y0)).^2 ./ (t/2).^2;
phi_r2= 1-((cos_theta).*(x-x0)+sin_theta.*(y-y0)).^2 ./ (l/2).^2;
phi_r= min(phi_r1,phi_r2);
phi_2D= max([phi_r,phi_c1,phi_c2],[],2) ;
phi_2D= gather(phi_2D);
end