
function [bins, ba] = eventsAlignedRasters2(st, eventTimes, trOrder, thisWindow, otherEvents, axRaster, axPSTH)

if isempty(trOrder)
    trOrder = 1:numel(eventTimes);
end

psthBinSize = 0.0001;
smoothWinStd = 0.01;
smoothWin = myGaussWin(smoothWinStd, 1/psthBinSize);

nTimes = length(eventTimes);
windowExp = thisWindow+smoothWinStd*5*[-1 1];
[ba, bins] = timestampsToBinned(st, eventTimes, psthBinSize, windowExp);
ba = ba(trOrder,:);

[tr,b] = find(ba);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything

thisPSTH = conv(nanmean(ba(~isnan(eventTimes),:))./psthBinSize, smoothWin, 'same');

axes(axRaster);
plot(rasterX, rasterY, 'k');
ylim([0 nTimes]);
xlim(thisWindow);
box off;
hold on;

% add other events nearby
reOrderEvent = eventTimes(trOrder);
for e = 1:length(otherEvents)
    otherEvent = otherEvents(e).times;
    x = WithinRanges(otherEvent, bsxfun(@plus, reOrderEvent, thisWindow), [1:nTimes]');
    [ii,trInds] = find(x);
    relTimes = otherEvent(ii)-reOrderEvent(trInds);
    plot(relTimes, trInds, otherEvents(e).icon, 'Color', otherEvents(e).color);
end

if ~isempty(axPSTH)
    axes(axPSTH);
    plot(bins, thisPSTH, 'k', 'LineWidth', 2.0); hold on;
    xlim(thisWindow);
    box off;
end