function finalShape = batch_union_all(barShapes, batchSize)
% batch_union_all
% Performs batch-wise union of polyshape objects using parallel processing,
% and recursively merges intermediate results into one final shape.
% A GUI progress bar is used for both batch and merge phases.
%
% Syntax:
%   finalShape = batch_union_all(barShapes, batchSize, mergeLevel)
%
% Inputs:
%   barShapes   - [1 x N] array of polyshape objects
%   batchSize   - Number of shapes per batch (default = 3)

%
% Output:
%   finalShape - Single polyshape after all unions
%
% Author: Zeyang Li

if nargin < 2, batchSize = 3; end
% if nargin < 3, mergeLevel = 20; end

finalShape=barShapes;
nBars = numel(barShapes);
iter = ceil(log(nBars)/log(batchSize));

p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= 4
    delete(gcp('nocreate'));
    parpool(4);
end

for i= 1:iter
nBars = numel(finalShape);
nBatch = ceil(nBars / batchSize);
groupUnion = repmat(polyshape(), 1, nBatch);

ss= sprintf('Batch Union:  (%d/%d)', i,iter);
WaitMsg = parfor_wait_ui(nBatch, 'Waitbar', true,'Title',ss,'ReportInterval', 5);
parfor k =1:nBatch
    iStart = (k - 1) * batchSize + 1;
   iEnd = min(k * batchSize, nBars);

  try
        groupUnion(k) = union(finalShape(iStart:iEnd));
    catch
        groupUnion(k) = polyshape();
  end

   WaitMsg.Send
   drawnow
   end

 finalShape =simplify(groupUnion);
 WaitMsg.Destroy;
end

end


% while  numel(polyList) > 1
% nBatch = ceil(nBars / batchSize);
% groupUnion = repmat(polyshape(), 1, nBatch);
% 
% 
% 
% ss= sprintf('Batch Union:  (/%d)',   nBatch) ;
% % WaitMsg = parfor_wait_ui(nBatch, 'Waitbar', true,'Title',ss,'ReportInterval', 10);
% 
% parfor k = 1:nBatch
%     iStart = (k - 1) * batchSize + 1;
%     iEnd = min(k * batchSize, nBars);
%     try
%         groupUnion(k) = union(barShapes(iStart:iEnd));
%     catch
%         groupUnion(k) = polyshape();
%     end
% 
%     % WaitMsg.Send
%     drawnow
% end
%   % WaitMsg.Destroy;
% 
% 
% % Proceed to final merge
% finalShape = recursive_merge(groupUnion, mergeLevel);
% end
