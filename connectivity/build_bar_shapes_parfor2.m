function barShapes = build_bar_shapes_parfor2(trussinfo)
% build_bar_shapes_parfor
% Constructs polyshape objects from xlist/ylist using parallel processing.
% Progress is shown via uiprogressdlg.
%
% Author: Zeyang Li
warning off;

% Ensure parallel pool
% p = gcp('nocreate');
% if isempty(p) || p.NumWorkers ~= 4
%     delete(gcp('nocreate'));
%     try;  parpool(4);
%     catch;    parpool(4);
%     end;
% end

[G, groupNames] = findgroups(trussinfo(:,5)');
groupCell = splitapply(@(idx){idx}, 1:length(trussinfo(:,5)), G);
n=length(groupCell );
barShapes = repmat(polyshape(), 1, n);
% Create GUI progress dialog
WaitMessage= parfor_wait_ui(n,  'Waitbar', true, 'ReportInterval', ceil(n/200),'Title', 'Generating polyshapes... ');
% Main loop
for i = 1:n
    id=groupCell{i};
    try
            th =  trussinfo( id(1),5)/2;
            coods =[ reshape( [trussinfo(id,1)' ; trussinfo(id,3)'; NaN(1, length(id))],1,[]) ;reshape( [trussinfo(id,2)' ; trussinfo(id,4)'; NaN(1, length(id))],1,[])];
            barShapes(i)=polybuffer(coods','lines', th);
    catch
         barShapes(i) = polyshape();
    end
  WaitMessage.Send;
  drawnow;
end
WaitMessage.Destroy;
end

