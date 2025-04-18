function opt=Zinterpo(opt,conkey)
[M,iK]= max(opt.V,[],2);
opt.Y=0;
for i= 1:opt.cand
se = (iK==i); 
opt.Y(se,i)=1; 
end

if conkey==2
opt.Y( M<0.35, [2,3])=1e-6;
opt.Y( M<0.35, [1])=1;
else
  opt.Y( M<0.35, :)= 1e-9; 
end
opt.Z= opt.Y;
opt.V= opt.Y;
end