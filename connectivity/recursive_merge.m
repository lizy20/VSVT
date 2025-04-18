function result = recursive_merge(polyList, level)
% recursive_merge
% Recursively merges a list of polyshape objects in groups of specified size (level).
% Merges are done in parallel per level and visualized using GUI progress.
%
% Syntax:
%   result = recursive_merge(polyList, level)
%
% Inputs:
%   polyList - Vector of polyshape objects
%   level    - Merge group size per round (e.g. 10)
%
% Output:
%   result   - Final polyshape after all recursive merges
%
% Author: Zeyang Li

if isempty(polyList)
    result = polyshape();
    return;
end

polyList = polyList(:);
roundNum = 0;

while numel(polyList) > 1
    roundNum = roundNum + 1;
    n = numel(polyList);
    nGroups = ceil(n / level);
    nextList = repmat(polyshape(), nGroups, 1);

 sf= sprintf('Merge Round:  (%d/%d)', ...
                roundNum, nGroups);
WaitMsg = parfor_wait_ui(nGroups, 'Waitbar', true,'Title',sf );
    parfor i = 1:nGroups
        iStart = (i - 1) * level + 1;
        iEnd = min(i * level, n);
        try
            nextList(i) = reducepolyunion(polyList(iStart:iEnd)) ;
        catch
            nextList(i) = polyshape();
        end
         WaitMsg.Send
    drawnow
    end
WaitMsg.Destroy;
    polyList = nextList;
end
result = polyList;

% 
% function result = recursive_merge(polyList, level)
% % recursive_merge
% % Recursively merges polyshapes in groups of 'level' using parallel union.
% % Shows GUI progress for each round.
% %
% % Author: Zeyang Li
% 
% if isempty(polyList)
%     result = polyshape();
%     return;
% end
% 
% polyList = polyList(:);
% roundNum = 0;
% 
% f = uifigure('Name', 'Recursive Merge');
% d = uiprogressdlg(f, ...
%     'Title', 'Recursive merging...', ...
%     'Message', 'Initializing...', ...
%     'Indeterminate', 'off');
% 
%  numel(polyList) > 1
%     roundNum = roundNum + 1;
%     n = numel(polyList);
%     nGroups = ceil(n / level);
%     nextList = repmat(polyshape(), nGroups, 1);
% 
%     dq2 = parallel.pool.DataQueue;
%     progress = 0;
%     total = nGroups;
% 
%     afterEach(dq2, @(~) updateProgress());
%      function updateProgress()
%         progress = progress + 1;
%         if mod(progress, 5) == 0 || progress == total
%             d.Value = progress / total;
%             d.Message = sprintf('Merge Round %d: %.1f%% (%d/%d)', ...
%                 roundNum, 100 * progress / total, progress, total);
%             drawnow;
%         end
%     end
% 
%     parfor i = 1:nGroups
%         iStart = (i - 1) * level + 1;
%         iEnd = min(i * level, n);
%         try
%             nextList(i) = union(polyList(iStart:iEnd));
%         catch
%             nextList(i) = polyshape();
%         end
%         send(dq2, []);
%         drawnow;
%     end
% 
%     polyList = nextList;
% end
% close(d); close(f);
% result = polyList;
% end


