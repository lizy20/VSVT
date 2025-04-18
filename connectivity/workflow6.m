function workflow6(filename,app)
%========Step 5:  Bound simplify=========================================
load([pwd, '\' , filename , '\geo_data5.mat']);

polyout = simplify(polyout);
border_thickness = 0.02*avglen; % 0 == No Border// 0.02-0.05recommand
simplifyTol =0;% 0 == No simple //  0.02 recommand

poly_total_bound = polyshape(nbounds(:,1), nbounds(:,2));
opt.v3 = sum(opt.V,2).* opt.xPhys(:);
Real_ele = opt.enidElem(opt.v3>0.3,:);
[MainShape, innerBorder] = build_region_and_border(opt.Node, Real_ele ,...
    border_thickness, simplifyTol);

polyout2 = union(polyout, innerBorder);
polyout_f = intersect(polyout2, MainShape);
polyout_f=simplify(polyout_f);
% ========Step 6: Plot=================================================
figure(); hold on; axis equal;
plot(polyout_f, 'FaceColor', 'black', 'FaceAlpha', 1, 'EdgeColor', 'none');
plot(nbounds(:,1), nbounds(:,2), 'r--', 'LineWidth', 1.5);
set(gca, 'YDir', 'normal');
title('Porous Structure');
%=================================================================
% Step 7:  volume constraint
area_out = polyout_f.area();
area_total = sum(opt.ElemArea(:));
fprintf(' Structure Practical Volume fraction: %.4f%%\n', area_out/area_total*100 );
myPrint(app, ' Structure Practical Volume fraction: %.4f%%\n', area_out/area_total*100 );
% Step 8:
var = setdiff(who, {'app'});
save([pwd, '\' , filename , '\geo_data6.mat'], var{:});
vol_diff = area_out/area_total - vol_lim;

end