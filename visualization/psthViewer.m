function psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
% psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
%
% GUI to display psths

params.smoothSize = 15; % in msec, stdev of gaussian smoothing filter
params.clusterIndex = 1;
params.rasterScale = 6; % height of ticks
params.window = window;

myData.spikeTimes = spikeTimes;
myData.clu = clu;
myData.eventTimes = eventTimes(:);
myData.trGroups = trGroups(:);
myData.params = params;
myData.clusterIDs = unique(clu);

f = figure; 

set(f, 'UserData', myData);
set(f, 'KeyPressFcn', @(f,k)psthViewerCallback(f, k));

psthViewerPlot(f)
end

function psthViewerPlot(f)

myData = get(f,'UserData');

% pick the right spikes
st = myData.spikeTimes(myData.clu==myData.clusterIDs(myData.params.clusterIndex));

% compute everything
[psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, myData.eventTimes, myData.params.window, 0.001);

% smooth the PSTH
gw = gausswin(round(myData.params.smoothSize*6),3);
smWin = gw./sum(gw);
psthSm = conv(psth, smWin, 'same');

% scale the raster ticks
rasterY(2:3:end) = rasterY(2:3:end)+myData.params.rasterScale;

% compute the tuning curve
trGroupLabels = unique(myData.trGroups);
nGroups = length(trGroupLabels);
tuningCurve = zeros(nGroups,2);
for g = 1:nGroups
    theseCounts = spikeCounts(myData.trGroups==trGroupLabels(g));
    tuningCurve(g,1) = mean(theseCounts);
    tuningCurve(g,2) = std(theseCounts)./sqrt(length(theseCounts));
end

 
% Make plots
subplot(3,1,1);
plot(bins, psthSm);
xlim(myData.params.window);
title(['cluster ' num2str(myData.clusterIDs(myData.params.clusterIndex))]);
xlabel('time (sec)');
ylabel('firing rate (Hz)');
makepretty;

subplot(3,1,2);
plot(rasterX,rasterY, 'k');
xlim(myData.params.window);
ylim([0 length(myData.eventTimes)+1]);
ylabel('event number');
xlabel('time (sec)');
makepretty;

subplot(3,1,3);
errorbar(trGroupLabels, tuningCurve(:,1), tuningCurve(:,2), 'o-');
xlabel('grouping variable value');
ylabel('average firing rate (Hz)');
makepretty;

end

function psthViewerCallback(f, keydata)

myData = get(f, 'UserData');

switch keydata.Key
    case 'rightarrow' % increment cluster index
        
        myData.params.clusterIndex = myData.params.clusterIndex+1;
        if myData.params.clusterIndex>length(myData.clusterIDs)
            myData.params.clusterIndex=1;
        end        
        
    case 'leftarrow' % decrement cluster index
        
        myData.params.clusterIndex = myData.params.clusterIndex-1;
        if myData.params.clusterIndex<1
            myData.params.clusterIndex=length(myData.clusterIDs);
        end
        
    case 'uparrow' % increase smoothing
        myData.params.smoothSize = myData.params.smoothSize*1.2;
        
    case 'downarrow' % decrease smoothing
        myData.params.smoothSize = myData.params.smoothSize/1.2;
        
end

set(f, 'UserData', myData);

% plot with new settings
psthViewerPlot(f)

end
