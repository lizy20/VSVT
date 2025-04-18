function barShapes = build_bar_shapes_parfor(xlist, ylist)
% build_bar_shapes_parfor
% Constructs polyshape objects from xlist/ylist using parallel processing.
% Progress is shown via uiprogressdlg.
%
% Author: Zeyang Li
warning off;
n = size(xlist, 1);
total = n;
barShapes = repmat(polyshape(), 1, n);

% Ensure parallel pool
p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= 4
    delete(gcp('nocreate'));
    try;  parpool(4);
    catch;    1;
    end;
end

% Create GUI progress dialog
WaitMessage= parfor_wait_ui(n,  'Waitbar', true, 'ReportInterval', ceil(n/200),'Title', 'Generating polyshapes... ');
% Main loop
parfor i = 1:n
    try
        barShapes(i) = polyshape(xlist(i,:), ylist(i,:));
    catch
        barShapes(i) = polyshape();
    end
  WaitMessage.Send;
  drawnow;
end
WaitMessage.Destroy;
end
