
function [avgPeriEventV, winSamps, periEventV, sortedLabels] = eventLockedAvg(V, t, eventTimes, eventLabels, calcWin)
% function [avgPeriEventV, periEventV] = eventLockedAvg(V, t, eventTimes, eventLabels, calcWin)
%
% Inputs: 
% - V: nTraces x nTimePoints
% - t: 1 x nTimePoints, the time points of each sample in V
% - eventTimes: 1 x nEvents, the time of each event. 
% - eventLabels: 1 x nEvents, the label of each event, e.g. the contrast
% value or some text label. If this is a cell array, the "tuning curve"
% will be plotted evenly spaced; if numeric array then these will be the
% x-axis values of the tuning curve
% - calcWin: 1 x 2, the start and end times relative to the event
%   - if calcWin is instead a vector with more than two elements, these are
%   taken to be the literal times around the event at which you want to
%   sample. E.g. you can use [-0.5 1.0] for the first version, or
%   [-0.5:0.1:1.0] for the second version (which allows you to specify the
%   sampling density)
%
% Outputs:
% - avgPeriEventV: nEventTypes x nTraces x nTimePoints, average temporal
% components across all events of each type
% - winSamps: labels for the time axis, relative to the event times
% - periEventV: nEvents x nTraces x nTimePoints, the temporal components around
% each event
% - sortedLabels: the labels of the rows of periEventV

t = t(:)'; % make row
eventTimes = eventTimes(:)';
[eventTimes, ii] = sort(eventTimes); %event times need to be in order for the interp function
sortedLabels = eventLabels(ii); % but that means the labels are rearranged too, and the user needs to know

nTraces = size(V,1);

eLabels = unique(eventLabels);
nConditions = length(eLabels);

if numel(calcWin)==2
    Fs = 1/median(diff(t));
    winSamps = calcWin(1):1/Fs:calcWin(2);
else
    winSamps = calcWin(:)';
end

periEventTimes = bsxfun(@plus, eventTimes', winSamps); % rows of absolute time points around each event
periEventV = zeros(nTraces, length(eventTimes), length(winSamps));
for s = 1:nTraces
    periEventV(s,:,:) = interp1(t, V(s,:), periEventTimes);
end

avgPeriEventV = zeros(nConditions, nTraces, length(winSamps));
for c = 1:nConditions
    if iscell(eventLabels)
        thisCondEvents = cellfun(@(x)strcmp(x,eLabels(c)),sortedLabels);
    else
        thisCondEvents = sortedLabels==eLabels(c);
    end
    avgPeriEventV(c,:,:) = squeeze(nanmean(periEventV(:,thisCondEvents,:),2));
end

periEventV = permute(periEventV, [2 1 3]);