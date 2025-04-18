function pointstransfer()


TR0=[1,0,-xl;
    0,1,-yl;
    0,0,1];
TR1=[1,0,xl;
    0,1,yl;
    0,0,1];
CR0=[cosd(ang), sind(ang),0; 
          -sind(ang), cosd(ang),0;
          0,0,1];
CR1 = [cosd(-ang), sind(-ang),0; 
          -sind(-ang), cosd(-ang),0;
          0,0,1];
      
SR0=[1/aniso,0,0; 0,1,0;0,0,1; ];
SR1=[aniso,0,0; 0,1,0;0,0,1; ];
Bs=SR0*CR0*TR0*boundary2;
boundary3=Bs([1,2],:)'; 

kp=[points, ones(length(points),1)]';
kps=SR0*CR0*TR0*kp;
kps2= kps([1,2],:)';
eps=0.005;

axis equal
vornum=length(C2);
vorbar=-1*(1000)*ones(vornum,20);
vorbardrawingrank=[0];
c=0;
vorthickness=[];
vormidpoints=[];
count2=0;
count1=0;
for i =1:vornum 
    vertex=V2([C2{i}],:);
    detect= isinf(vertex);
    if any(detect)==1
        vertex=vertex(:);
        vertex(find(detect))=1e5;
        vertex=reshape(vertex, [],2 );
    end
    IDX = convhull(vertex);
    vertex= vertex(IDX,1:2);
    vertex=vertex';
    vorbar(i,[1:length(vertex(:))])=vertex(:);
    vorbarnum(i,1)=length(IDX);
    count1=count1+length(IDX)-1;
    if count1>500
        count1=0;
        vorbardrawingrank=[vorbardrawingrank,i];
    end
    for j= 1:length(IDX)-1%-1
        count2=count2+1;
        vormidpoints=[ vormidpoints;[(vorbar(i,2*j-1)+vorbar(i,2*j+1))/2,(vorbar(i,2*j)+vorbar(i,2*j+2))/2]];
        vorbarinfo(count2,:) = (vorbar(i,[2*j-1 : 2*j+2]));
    end
end

cc=[vorbarinfo(:,1) , vorbarinfo(:,2), ones(length(vorbarinfo),1) ]';
cc = (TR1*CR1*SR1*cc)';
vorbarinfo(:,1)=cc(:,1); vorbarinfo(:,2)= cc(:,2);

cc=[vorbarinfo(:,3), vorbarinfo(:,4), ones(length(vorbarinfo),1)]';
cc = (TR1*CR1*SR1*cc)';
vorbarinfo(:,3)= cc(:,1); vorbarinfo(:,4)=cc(:,2);
time=toc;
ob=ismembertol(vorbarinfo(:,1),boundary(:,1),0.001)+ ismembertol(vorbarinfo(:,3),boundary(:,1),0.001)+...
   ismembertol(vorbarinfo(:,2),boundary(:,2),0.001)+ ismembertol(vorbarinfo(:,4),boundary(:,2),0.001);
vorbarinfo(ob>=2,:)=[]; %1
disp(['泰勒多边形用时' num2str(time)] );
end