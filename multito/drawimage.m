
function []= drawimage(opt,app)

Vertices=opt.Node;
if size(opt.V,2)>1
[Color,~] = ConvertColors(opt.V);
I =reshape(Color, [],3);
cla(app.UIAxes2);
hold(app.UIAxes2, 'on'); 
title(app.UIAxes2, 'Anisotropy');
c1 = plot(app.UIAxes2, nan, nan, 's', 'MarkerSize', 10, ...
    'MarkerFaceColor', [   0.5647    0.7961    0.9843], 'MarkerEdgeColor', 'k', 'DisplayName', ['Aiso=', num2str(opt.ml(1))]);
c2 = plot(app.UIAxes2, nan, nan, 's', 'MarkerSize', 10, ...
    'MarkerFaceColor', [ 0.9647    0.7922    0.8980], 'MarkerEdgeColor', 'k', 'DisplayName', ['Aiso=', num2str(opt.ml(2))]);
c3 = plot(app.UIAxes2, nan, nan, 's', 'MarkerSize', 10, ...
    'MarkerFaceColor', [  0.5882    0.8000    0.7961], 'MarkerEdgeColor', 'k', 'DisplayName', ['Aiso=', num2str(opt.ml(3))]);
lgd = legend(app.UIAxes2,"show" );    

patch(app.UIAxes2, 'Faces',padcellrows(opt.Element), 'Vertices',Vertices,'FaceVertexCData',I ,...
     'Facecolor', 'flat','EdgeColor','none','HandleVisibility', 'off');
axis(app.UIAxes2, 'equal'); 
drawnow 
hold(app.UIAxes2, 'off');

else
cla(app.UIAxes2);
 patch(app.UIAxes2,'Faces',padcellrows(opt.Element),'Vertices',Vertices,'FaceVertexCData',opt.V,...
     'Facecolor', 'flat','EdgeColor','none');
axis(app.UIAxes2, 'equal');  drawnow
end

patch(app.UIAxes3,'Faces',padcellrows(opt.Element),'Vertices',Vertices,'FaceVertexCData',opt.xPhys,...
     'Facecolor', 'flat' ,'EdgeColor','none'); 
clim([0 1]);colorbar(app.UIAxes3);  axis(app.UIAxes3, 'equal');  drawnow

patch(app.UIAxes4,'Faces',padcellrows(opt.Element),'Vertices',Vertices,'FaceVertexCData',opt.alpha,...
     'Facecolor', 'flat','EdgeColor','none');
clim([-90 90]); colorbar(app.UIAxes4);  axis(app.UIAxes4, 'equal');  drawnow

patch(app.UIAxes5, 'Faces',padcellrows(opt.Element),'Vertices',Vertices,'FaceVertexCData',opt.vx2,...
     'Facecolor', 'flat'); 
clim(app.UIAxes5, [0,1]); 
axis(app.UIAxes5, 'equal');  drawnow



end