


function [rocOut, tpr, fpr] = rocArea(values1, values2)
% function [rocOut, tpr, fpr] = rocArea(values1, values2)
%
% "values2" should have the larger mean, if you want it to return a value
% greater than 0.5. 

if isempty(values1) && isempty(values2)
    rocOut = NaN;
    tpr = [];
    fpr = [];
    return;
end

% check that inputs are row vectors
if size(values1,1)>size(values1,2)
    values1 = values1';
end
if size(values2,1)>size(values2,2)
    values2 = values2';
end

[tpr, fpr, th] = roc([zeros(1,length(values1)) ones(1,length(values2))], [values1 values2]);
tpr = [tpr 1];
fpr = [fpr 1];

% old definition of area: only correct in the limit of infinite samples!
% Otherwise overestimates area.
% rocOut = sum(diff(fpr).*tpr(2:end));

% new definition of area: exactly calculates, so that AUROC is commutative.
dfpr = diff(fpr);
dtpr = diff(tpr);
% the first term is the little triangles under each line segment connecting
% adjacent points. the second term is the rectangles under these triangles.
rocOut = 0.5*sum(dfpr.*dtpr) + sum(tpr(2:end-1).*dfpr(2:end));

% if mean(values1)<mean(values2)
%     rocOut = sum(diff(fpr).*tpr(2:end));
% else
%     rocOut = 1-sum(diff(fpr).*tpr(2:end));
% end