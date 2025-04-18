function [De_t, dDe_t] = DH_fit2D(E0,Emin,x,m,DNN)
De_t=zeros(9,length(x),m);
dDe_t=zeros(9,length(x),m);
nele=length(x);
u=0.3;
E0=1*E0;
De=E0/(1-u^2)*[1,u,0;u,1,0;0,0,(1-u)/2];
Demin=1e-6/(1-u^2)*[1,u,0;u,1,0;0,0,(1-u)/2];

for aid=1:m          
n1=DNN{aid,1}.net;
n2=DNN{aid,4}.net;
n3=DNN{aid,2}.net; 
n4=DNN{aid,3}.net;

D11 =   Demin(1,1)+(E0 )*n1(x')' ;
D12 =   Demin(1,2)+( E0 )*n2(x')'  ;        
D21=  D12;     
D22 =    Demin(2,2)+ ( E0 )*n3(x')'  ; 
D33 =  Demin(3,3) + ( E0 )*n4(x')' ; 
nele=length(x);


D11(x>=1)=De(1,1); D12(x>=1)=De(1,2); D21(x>=1)= De(1,2);
D22(x>=1)=De(2,2); D33(x>=1)=De(3,3); 
De_t(:,:,aid) = [D11,D21,zeros(nele,1), D12,D22,zeros(nele,1),... 
    zeros(nele,1),zeros(nele,1),D33]';

d_D11 = ( E0 )*calgradient(n1,x);
d_D21 =( E0  ) *calgradient(n2,x);
d_D12 =  d_D21;
d_D22 = ( E0  )*calgradient(n3,x);
d_D33 = ( E0 )*calgradient(n4,x);
d_D11(x>=1)=0; d_D12(x>=1)=0; d_D21(x>=1)= 0;
d_D22(x>=1)=0; d_D33(x>=1)=0; 

dDe_t(:,:,aid) = [d_D11,d_D21,zeros(nele,1), d_D12,d_D22,zeros(nele,1),... 
    zeros(nele,1),zeros(nele,1),d_D33]';

end


end