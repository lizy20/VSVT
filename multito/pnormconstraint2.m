
function [v3, dv3dV]=pnormconstraint2(V)
pq= -16 ; %-10 
[n,m]=size(V);
vec=sum(V,2);
v3= - norm(vec(:), pq)+0.99; %0.95
dv3dV=-1*repmat( ((sum(vec(:).^pq))^(1/pq-1)) * (vec(:).^(pq-1)), [1,m]);
end