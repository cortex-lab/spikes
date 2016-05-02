function [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(spikeTimes, eventTimes, window, psthBinSize)
% function [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(spikeTimes, eventTimes, window)
%
% Fast computation of psth and spike counts in a window relative to some
% events. Also returns rasters you can plot. 
%
% Notes on inputs:
% - eventTimes is nEvents x 1
% - window is length 2, e.g. [-0.1 0.3] for a window from -0.1 to +0.3 sec
% relative to events. 
%
% Notes on outputs:
% - psth can be plotted with plot(bins, psth);
% - rasters can be plotted with plot(rasterX, rasterY);
% - spikeCounts is nEvents x 1, where each entry is the number of spikes
% that occurred within the window around that event. 


nRanges = size(eventTimes,1);
ranges = [eventTimes+window(1) eventTimes+window(2)];
rangeLabel = 1:nRanges;

% wr will have one entry per spike, corresponding to an integer identifying
% which range the spike was in
wr = WithinRanges(spikeTimes, ranges, rangeLabel', 'vector');

stIn = spikeTimes(wr>0); wr = wr(wr>0); % pick just the spikes and range indices that actually were within a range
stRelToEvent = stIn-eventTimes(wr); % subtract the event time corresponding to each spike


[psth,bins] = hist(stRelToEvent, [window(1):psthBinSize:window(2)]);


[rasterX,yy] = rasterize(stRelToEvent); 
rasterY = yy+reshape(repmat(wr,3,1),1,length(wr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything


spikeCounts = zeros(1, nRanges);

% inner diff gives the indices when you go from spikes of one range to
% spikes of the next. The distance between these breaks are the spike counts
% Find gives you the indices of these breaks, outer diff computes how many in
% each.
spikeCounts(ismember(1:nRanges,wr)) = diff(find(diff([0 wr nRanges+1])>0)); 