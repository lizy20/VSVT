function [v2, dv2dV]=pnormconstraint(V)
pq = 12;
[n,m]=size(V);
vec=sum(V,2);
v2= norm(vec(:), pq) - 1 ;
dv2dV=repmat( ((sum(vec(:).^pq))^(1/pq-1)) * (vec(:).^(pq-1)), [1,m]);
end