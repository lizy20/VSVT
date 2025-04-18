
function run_geo_main(L, W,  Nx, Ny, VolFrac, PL, PU,app)

startBusyLamp(app); 
pause(0.05);
M = [ L,W, Nx,Ny, VolFrac, PL, PU];
vf=M(5);

se=app.BoundaryConditionDropDown.Value;
name= app.BoundaryConditionDropDown.Items;
filename= [name{se}, '_' ,num2str(M(6)),'-',num2str(M(7)), ...
                '_', num2str(M(3)), '-',num2str(M(4)),  '_', num2str(vf) ];
% filename='Cantilever_0.3-0.8_50-30_0.5';
CStatus = app.getStepCheckVector();
if CStatus(1)==1; workflow1(filename,app); end
if CStatus(2)==1;workflow2(filename,app); end
if CStatus(3)==1;workflow3(filename,app);end
if CStatus(4)==1;workflow4(filename,app);end
if CStatus(5)==1;workflow5(filename,app);end
if CStatus(6)==1;workflow6(filename,app);end
stopBusyLamp(app); 
end
