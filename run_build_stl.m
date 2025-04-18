function run_build_stl(L, W,  Nx, Ny, VolFrac, PL, PU,app)
M = [ L,W, Nx,Ny, VolFrac, PL, PU];
vf=M(5);
startBusyLamp(app); 
se=app.BoundaryConditionDropDown.Value;
name= app.BoundaryConditionDropDown.Items;
filename= [name{se}, '_' ,num2str(M(6)),'-',num2str(M(7)), ...
                '_', num2str(M(3)), '-',num2str(M(4)),  '_', num2str(vf) ];
load([pwd, '\' , filename , '\geo_data6.mat']);
polyout_f=simplify(polyout_f);
tr =triangulation(polyout_f);
stlwrite(tr,[pwd,'\' ,filename,'\', filename ,'.stl'], 'binary'); % polyfunc toolbox built-in func
stopBusyLamp(app); 
end