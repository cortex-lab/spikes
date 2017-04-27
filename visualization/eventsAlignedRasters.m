

function eventsAlignedRasters(spikeTimes, eventTimes, eventNames, windows)
% function eventsAlignedRasters(eventTimes, spikeTimes, eventNames, windows)
%
% plots rasters of the spikes aligned to each of multiple events, with
% relative times of the other events indicated.
% - spikeTimes is a vector
% - other three arguments are cell arrays, one for each event. Assumes that
% all three have the same number of elements.

psthBinSize = 0.0001;
nEv = length(eventTimes);
nTimes = length(eventTimes{1});
smoothWinStd = 0.01;
smoothWin = myGaussWin(smoothWinStd, 1/psthBinSize);

for e = 1:nEv
    subplot(3,nEv,[e e+nEv]); hold off
        
    for e2 = 1:nEv
        plot(eventTimes{e2}-eventTimes{e}, 1:nTimes, 'o'); hold on;
    end
    
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
    plot(bins, conv(nanmean(ba(~isnan(eventTimes{e}),:))./psthBinSize, smoothWin, 'same'), 'k', 'LineWidth', 2.0); hold on;
    xlim(windows{e});
    box off;
    if e==1
        ylabel('firing rate (sp/s)');
    end
        
    xlabel(sprintf('time from %s (s)', eventNames{e}));
end
equalizeAxes(axFR,'y');
co = get(gca, 'ColorOrder');    
yl = ylim();
for e = 1:nEv
    axes(axFR(e));    
    plot([0 0], yl, '--', 'Color', co(e,:), 'LineWidth', 2.0);
    ylim(yl);
end