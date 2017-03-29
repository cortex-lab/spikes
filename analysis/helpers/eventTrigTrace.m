

function [mn, stdErr, winSamps, allTraces] = eventTrigTrace(traceT, trace, eventTimes, win, varargin)
% function [mn, stdev, allTraces] = eventTrigTrace(traceT, trace, eventTimes, win[, ops])
%
% ops.Fs allows you to pick a different timescale of the average than of
% the original trace, otherwise will use from the original trace.
%
% ops.doSubBsl, boolean, will subtract from each trace the mean of the
% before time-zero data for each, before averaging
ops = [];
if ~isempty(varargin)
    ops = varargin{1};
end
ops.doSubBsl = pick(ops, 'doSubBsl', 'def', false);
ops.Fs = pick(ops, 'Fs', 'def', 1/mean(diff(traceT)));

nEvents = length(eventTimes);
        
winSamps = win(1):1/ops.Fs:win(2);
periEventTimes = bsxfun(@plus, eventTimes, winSamps); % rows of absolute time points around each event

allTraces = interp1(traceT, trace, periEventTimes);

if ops.doSubBsl
    bslMn = nanmean(allTraces(:,winSamps<0),2);
    allTraces = bsxfun(@minus, allTraces, bslMn);
end

mn = nanmean(allTraces);
stdErr = std(allTraces,0,1,'omitnan')./sqrt(nEvents);