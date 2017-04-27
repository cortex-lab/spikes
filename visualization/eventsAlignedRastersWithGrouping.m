

function eventsAlignedRastersWithGrouping(spikeTimes, eventTimes, params)
% function eventsAlignedRasters(eventTimes, spikeTimes, params)
%
% plots rasters of the spikes aligned to each of multiple events, with
% relative times of the other events indicated.
% - spikeTimes is a vector
% - eventTimes is cell arrays, one for each event. 
% - params has:
%   - eventNames, cell array with names for events
%   - windows, cell array with start and end time around each event to use
%   - eventGroups, cell array with a vector length nEvents in each cell,
%   the unique values of which specify the different groups to split by
%   - colors, cell array, one for each event type, contains the colors to
%   use for each unique value
%
% Will sort each plot independently, based on the groupings for that event
%
% In this version, can have a different number of events in each plot,
% doesn't matter 

psthBinSize = 0.0001;
nEv = length(eventTimes);
smoothWinStd = 0.01;
smoothWin = myGaussWin(smoothWinStd, 1/psthBinSize);

windows = params.windows;
eventNames = params.eventNames;

for e = 1:nEv
    
    nTimes = length(eventTimes{e});
    
    uGroups = unique(params.eventGroups{e});
    nG = numel(uGroups);
    [groupLabels, ii] = sort(params.eventGroups{e});
    
    
    subplot(3,nEv,[e e+nEv]); hold off
        
    for g = 1:nG
        theseTr = find(groupLabels==uGroups(g) & ~isnan(eventTimes{e}));
        plot(zeros(size(theseTr)), theseTr, 'o', 'Color', params.colors{e}(g,:));
        hold on;
    end
    
    % TODO: take care of plotting the other events relative to this one
%     for e2 = 1:nEv
%         plot(eventTimes{e2}-eventTimes{e}, 1:nTimes, 'o'); hold on;
%     end
    
    windowExp = windows{e}+smoothWinStd*5*[-1 1];
    [ba, bins] = timestampsToBinned(spikeTimes, eventTimes{e}, psthBinSize, windowExp);        
    
    [tr,b] = find(ba);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
    
    plot(rasterX, rasterY, 'k');  
    ylim([0 nTimes]);
    
    if e==1
        ylabel('event number');
    end
    xlim(windows{e});
    box off;
    
    axFR(e) = subplot(3,nEv,e+2*nEv); hold off;
    
    for g = 1:nG
        theseTr = groupLabels==uGroups(g) & ~isnan(eventTimes{e});
    
        plot(bins, conv(nanmean(ba(theseTr,:))./psthBinSize, smoothWin, 'same'), ...
            'Color', params.colors{e}(g,:), 'LineWidth', 2.0); 
        hold on;
    end
    xlim(windows{e});
    box off;
    if e==1
        ylabel('firing rate (sp/s)');
    end
        
    xlabel(sprintf('time from %s (s)', eventNames{e}));
end
equalizeAxes(axFR,'y');
yl = ylim();
for e = 1:nEv
    axes(axFR(e));    
    plot([0 0], yl, 'k--', 'LineWidth', 2.0);
    ylim(yl);
end