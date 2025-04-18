function [handle,map] = InitialPlot(opt)
Tri = zeros(length([opt.Element{:}])-2*opt.nele,3);
map = zeros(size(Tri,1),1); index=0;
for el = 1:opt.nele
  for enode = 1:length(opt.Element{el})-2
    map(index+1) = el;
    Tri(index+1,:) = opt.Element{el}([1,enode+1,enode+2]);
    index = index + 1;
  end
end
I=[zeros(opt.nele,3),zeros(opt.nele,1)];
Vertices= opt.Node;
handle = patch(gca,'Faces',Tri,'Vertices',opt.Node,'FaceVertexCData',...
               I(map,1:3),'FaceColor','flat','EdgeColor','none');
% axis equal; axis off;
end