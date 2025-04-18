function opt = updateinfo(opt)
% opt.Z = opt.V;
% opt.Y= opt.V 
% [opt.W, dwdy,opt.V,dVdy] = MatIntFnc(opt.Y, 0.5, 1,'SIMP-H',3);
opt.xPhys = opt.P*opt.x(:);  
opt.vx2 = sum(opt.V,2).*opt.xPhys(:);
disp(  [ 'final volume:' , num2str( sum(opt.vx2(:))/sum(opt.ElemArea))] );
end