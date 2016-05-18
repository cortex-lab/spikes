

function psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
% function psthViewer(spikeTimes, clu, eventTimes, window, trGroups)
%
% Controls:
% - c: dialog box to pick a new cluster ID number
% - t: toggle showing psth traces for each grouping variable or just the
% overall. If showing just overall, raster is sorted chronologically. If
% showing by grouping variable, raster is sorted by that variable.
% - r: select a new range within which to count spikes for the tuning curve
% -
%
% TODO:
% - if plotting all traces, color raster and sort
% - don't replot traces just change X,Y data for speed and for keeping zoom
% levels
% - indicate on tuning curve which color is which
% - add ability to switch from graded to distinguishable color scheme

params.smoothSize = 15; % in msec, stdev of gaussian smoothing filter
params.clusterIndex = 1;
params.rasterScale = 6; % height of ticks
params.window = window;
params.showAllTraces = false;
params.showErrorShading = false;
params.startRange = window(1);
params.stopRange = window(2);
params.binSize = 0.001;

myData.spikeTimes = spikeTimes;
myData.clu = clu;
myData.eventTimes = eventTimes(:);
myData.trGroups = trGroups(:);
myData.clusterIDs = unique(clu);
myData.trGroupLabels = unique(myData.trGroups);
myData.nGroups = length(myData.trGroupLabels);

params.colors = copper(myData.nGroups); params.colors = params.colors(:, [3 2 1]);

myData.params = params;

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
%[psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, myData.eventTimes, myData.params.window, 0.001);
[psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(st, myData.eventTimes, myData.params.window, myData.params.binSize);
trGroupLabels = myData.trGroupLabels;
nGroups = myData.nGroups;
inclRange = bins>myData.params.startRange & bins<=myData.params.stopRange;
spikeCounts = sum(ba(:,inclRange),2)./(myData.params.stopRange-myData.params.startRange);


% PSTH smoothing filter
gw = gausswin(round(myData.params.smoothSize*6),3);
smWin = gw./sum(gw);

% smooth ba
baSm = conv2(smWin,1,ba', 'same')'./myData.params.binSize;

% construct psth(s) and smooth it
if myData.params.showAllTraces
    psthSm = zeros(nGroups, numel(bins));
    if myData.params.showErrorShading
        stderr = zeros(nGroups, numel(bins));
    end
    for g = 1:nGroups
        psthSm(g,:) = mean(baSm(myData.trGroups==trGroupLabels(g),:));
        if myData.params.showErrorShading
            stderr(g,:) = std(baSm)./sqrt(size(baSm,1));
        end
    end
else
    
    psthSm = mean(baSm);
    if myData.params.showErrorShading
        stderr = std(baSm)./sqrt(size(baSm,1));
    end
    
end

% compute raster
if myData.params.showAllTraces
    [~, inds] = sort(myData.trGroups);
    [tr,b] = find(ba(inds,:));
else
    [tr,b] = find(ba);
end
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything

% scale the raster ticks
rasterY(2:3:end) = rasterY(2:3:end)+myData.params.rasterScale;

% compute the tuning curve
tuningCurve = zeros(nGroups,2);
for g = 1:nGroups
    theseCounts = spikeCounts(myData.trGroups==trGroupLabels(g));
    tuningCurve(g,1) = mean(theseCounts);
    tuningCurve(g,2) = std(theseCounts)./sqrt(length(theseCounts));
end


% Make plots
colors = myData.params.colors;
subplot(3,1,1); hold off;
if myData.params.showAllTraces
    for g = 1:nGroups
        plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
        hold on;
    end
else
    plot(bins, psthSm);
end
xlim(myData.params.window);
title(['cluster ' num2str(myData.clusterIDs(myData.params.clusterIndex))]);
xlabel('time (sec)');
ylabel('firing rate (Hz)');
yl = ylim();
hold on;
plot(myData.params.startRange*[1 1], yl, 'k--');
plot(myData.params.stopRange*[1 1], yl, 'k--');
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
        
    case 'e' % whether to show standard error as shading
        myData.params.showErrorShading = ~myData.params.showErrorShading;
        
    case 't' % whether to plot the psth trace for each condition or just the overall one
        myData.params.showAllTraces = ~myData.params.showAllTraces;
        
    case 'r'
        ax = subplot(3,1,1); title('click start and stop of range')
        %         [startRange,~] = ginput();
        %         [stopRange,~] = ginput();
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        myData.params.startRange = q(1,1);
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        myData.params.stopRange = q(1,1);
        
    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            myData.params.clusterIndex = ind;
        end
        
end

set(f, 'UserData', myData);

% plot with new settings
psthViewerPlot(f)

end


function makepretty()
% set some graphical attributes of the current axis

set(get(gca, 'XLabel'), 'FontSize', 17);
set(get(gca, 'YLabel'), 'FontSize', 17);
set(gca, 'FontSize', 13);

set(get(gca, 'Title'), 'FontSize', 20);

ch = get(gca, 'Children');

for c = 1:length(ch)
    thisChild = ch(c);
    if strcmp('line', get(thisChild, 'Type'))
        if strcmp('.', get(thisChild, 'Marker'))
            set(thisChild, 'MarkerSize', 15);
        end
        if strcmp('-', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2.0);
        end
    end
end

end