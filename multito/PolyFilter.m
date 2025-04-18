%------------------------------ PolyTop ----------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyTop: A Matlab %
% implementation of a general topology optimization framework using       %
% unstructured polygonal finite element meshes", Struct Multidisc Optim,  %
% DOI 10.1007/s00158-011-0696-x                                           %
%-------------------------------------------------------------------------%
function [P] = PolyFilter(fem,R)
if R<0, P = speye(fem.nele); return; end %P is set to identity when R<0

% if fem.nonreg ==0
% [Den,Dens] = filtering2d(fem.nelx, fem.nely, fem.nele, R/fem.Elex);
% % P = Den./repmat(Dens,[1,fem.nele]); 
% 
% P = spdiags(1 ./ Dens, 0, fem.nele, fem.nele) * Den; 
% else

ElemCtrd = zeros(fem.nele,2);
for el = 1:fem.nele        %Compute the centroids of all the elements
  vx=fem.Node(fem.Element{el},1); vy=fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  A = 0.5*sum(temp);
  ElemCtrd(el,1) = 1/(6*A)*sum((vx+vx([2:end 1])).*temp);
  ElemCtrd(el,2) = 1/(6*A)*sum((vy+vy([2:end 1])).*temp);
end
if ~isfield(fem,'SElem'), SElem = [];
else SElem = fem.SElem; end
% ElemCtrd = ElemCtrd + 1e-2 * randn(size(ElemCtrd)); 
%Perturbation to prevent non-convergence due to integer filtering radii

[d] = DistPntSets(ElemCtrd,ElemCtrd,R*1.00001,SElem);  %Obtain distance values & indices
% P = sparse(d(:,1),d(:,2),1-d(:,3)/R);    %Assemble the filtering matrix
% P = sparse(d(:,1),d(:,2),(1-d(:,3)/R).^2);  %polynomial filter  ， integer R may non-convergence

sigma = R / 3;    %  d = R ， w ≈ exp(-4.5) ≈ 0.011  gauss filter
idx = d(:,3) <= R;  % 
w = exp( - (d(idx,3).^2) / (2 * sigma^2) );
P = sparse(d(idx,1), d(idx,2), w);

P = spdiags(1./sum(P,2),0,fem.nele,fem.nele)*P;
% end
%---------------------------------- COMPUTE DISTANCE BETWEEN TWO POINT SETS

function [d] = DistPntSets(PS1,PS2,R,SElem)
d = cell(size(PS1,1),1);
for el = 1:size(PS1,1)       %Compute the distance information
  if ismember(el,SElem)
    d{el} = [el,el,0];
  else
    dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
    [I,J] = find(dist<=R);   %Find the indices for distances less that R
    [~,K] = setdiff(I,SElem);
    d{el} = [I(K),J(K)+(el-1),dist(I(K))];
  end
end
d = cell2mat(d);             %Matrix of indices and distance value
%-------------------------------------------------------------------------%