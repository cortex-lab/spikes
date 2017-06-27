
function plotMat = plotAsPhase3(data,xc, yc)
% data is a vector of length nChans

[data, xc, yc] = replaceMissingSitesP3(data, xc, yc);

data = data(:); % is column now

% turn xc into inds
uxc = unique(xc);
for u = 1:length(uxc)
    xc(xc==uxc(u))=u;
end

% turn yc into inds - here *assuming* a regular scale
% ycDiff = diff(yc); ycDiff = min(ycDiff(ycDiff>0));
ycDiff = 20;
yc = round(yc./ycDiff)+1;


plotMat = NaN(length(uxc)*2, max(yc)+1);

for d = 1:length(data)
    plotMat((xc(d)-1)*2+[1 2], (yc(d)-1)+[1 2]) = data(d);
end


h = imagesc((0:size(plotMat,1)-1)*ycDiff, (0:size(plotMat,2)-1)*ycDiff, plotMat');
ax = get(h, 'Parent');
% set(ax, 'XTick', [], 'YTick', []); 
set(ax, 'YDir', 'normal');
axis image
% truesize(gcf, [size(plotMat,1)*2 size(plotMat,2)*2]);
set(h,'alphadata',~isnan(plotMat)')
axis off
% colormap hot;



%--- old version:

% nCh = length(data);
% 
% plotMat = NaN(nCh,5);
% 
% thisCol = mod(1:nCh,4)==2;
% plotMat(1:4:end,1) = data(thisCol);
% plotMat(1:4:end,2) = data(thisCol);
% plotMat(2:4:end,1) = data(thisCol);
% plotMat(2:4:end,2) = data(thisCol);
% 
% thisCol = mod(1:nCh,4)==1;
% plotMat(1:4:end,3) = data(thisCol);
% plotMat(1:4:end,4) = data(thisCol);
% plotMat(2:4:end,3) = data(thisCol);
% plotMat(2:4:end,4) = data(thisCol);
% 
% thisCol = mod(1:nCh,4)==0;
% plotMat(3:4:end,2) = data(thisCol);
% plotMat(3:4:end,3) = data(thisCol);
% plotMat(4:4:end,2) = data(thisCol);
% plotMat(4:4:end,3) = data(thisCol);
% 
% thisCol = mod(1:nCh,4)==3;
% plotMat(3:4:end,4) = data(thisCol);
% plotMat(3:4:end,5) = data(thisCol);
% plotMat(4:4:end,4) = data(thisCol);
% plotMat(4:4:end,5) = data(thisCol);