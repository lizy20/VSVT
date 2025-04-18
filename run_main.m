% ====================================================================
% % % MIT License
% % % Copyright (c) 2025 [Zeyang Li,M2O group, School of engineering, Cardiff University]
%  For probable citaion for the codes,  please refer to 
%                                                 https://doi.org/10.1016/j.cma.2024.117378
%-----------------------------------------------------------------------------------------------------------------------
% % % Permission is hereby granted, free of charge, to any person obtaining a copy
% % % of this software and associated documentation files (the "Software"), to deal
% % % in the Software without restriction, including without limitation the rights
% % % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % % copies of the Software, and to permit persons to whom the Software is
% % % furnished to do so, subject to the following conditions:
% % % 
% % % The above copyright notice and this permission notice shall be included in all
% % % copies or substantial portions of the Software.
%===================================================================
function run_main(L, W,  Nx, Ny, VolFrac, PL, PU,rmin,Fkey,app)
M = [ L,W, Nx,Ny, VolFrac, PL, PU];
vf=M(5);
startBusyLamp(app); 
% M = [80, 80, 20,20, 0.3, 0.2,0.5];
% rmin=0.5;   Fkey=1 ;vf=0.5；
se=app.BoundaryConditionDropDown.Value;
name= app.BoundaryConditionDropDown.Items;
filename= [name{se}, '_' ,num2str(M(6)),'-',num2str(M(7)), ...
                '_', num2str(M(3)), '-',num2str(M(4)),  '_', num2str(vf) ];
switch (app.Switch.Value)
    case(1);  continous(filename, M, vf, 3, [1,1.5,3], Fkey, rmin,app);
    case(0);  continous(filename, M, vf, 1, [1], Fkey, rmin,app); end

postfunction(filename,app);
stopBusyLamp(app); 
end
% %==================================

