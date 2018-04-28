
function t = totalMinDiffs(b, predictTo, predictFrom)

myPred = [predictFrom(:) ones(size(predictFrom(:)))]*b;
md = findMinDiffs(myPred, predictTo).^2;
% throw away the worst N where N is the number of excess elements in the
% longer vector
mds = sort(md);
md = mds(1:min(length(predictTo), length(predictFrom)));
    
t = sum(md);

