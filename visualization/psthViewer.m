

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
% - add ability to switch from graded to distinguishable color scheme, also
% cyclical
% - add support for a second grouping variable - in that case, should use a
% different graded color scheme for each, overlay the traces in the tuning
% curve view, and also provide a 2-d image of tuning curve
% - add support for plot labels (traces in psth, x-axis in tuning curve)
% - option to change filter type
% XX expand calc window to avoid falling off at the edges
% - option for error bars
% - fix bug where the raster from the last unit stays when there are no new
% spikes to plot

fprintf(1, 'Controls:\n')
fprintf(1, '- left/right arrow: select previous/next cluster\n')
fprintf(1, '- up/down arrow: change smoothing of psth curves\n')
fprintf(1, '- c: dialog box to pick a new cluster ID number\n')
fprintf(1, '- t: toggle showing psth traces for each grouping variable or just the\n')
fprintf(1, 'overall. If showing just overall, raster is sorted chronologically. If\n')
fprintf(1, 'showing by grouping variable, raster is sorted by that variable.\n')
fprintf(1, '- r: select a new range within which to count spikes for the tuning curve\n')


params.smoothSize = 15; % in msec, stdev of gaussian smoothing filter
params.clusterIndex = 1;
params.rasterScale = floor(numel(eventTimes)/100); % height of ticks
params.window = window;
params.showAllTraces = false;
params.showErrorShading = false;
params.startRange = window(1);
params.stopRange = window(2);
params.binSize = 0.001;
params.filterType = 1; 
params.smWin = genSmWin(params);

myData.spikeTimes = spikeTimes;
myData.clu = clu;

if ~issorted(eventTimes)
    [eventTimes, ii] = sort(eventTimes);
    trGroups = trGroups(ii);
end

myData.eventTimes = eventTimes(:);
myData.trGroups = trGroups(:);
myData.clusterIDs = unique(clu);
myData.trGroupLabels = unique(myData.trGroups);
myData.nGroups = length(myData.trGroupLabels);
myData.plotAxes = [];

params.colors = copper(myData.nGroups); params.colors = params.colors(:, [3 2 1]);

myData.params = params;

f = figure; f.Color = 'w';
set(f, 'Renderer', 'painters')

set(f, 'UserData', myData);
set(f, 'KeyPressFcn', @(f,k)psthViewerCallback(f, k));

psthViewerPlot(f)
end

function psthViewerPlot(f)
% fprintf(1,'plot with fig %d\n', get(f,'Number'));
myData = get(f,'UserData');
p = myData.params;

% pick the right spikes
st = myData.spikeTimes(myData.clu==myData.clusterIDs(myData.params.clusterIndex));

plotWindow = p.window;
calcWindow = plotWindow+p.smoothSize*3/1000*[-1 1];

% compute everything
%[psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st, myData.eventTimes, myData.params.window, 0.001);
[psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(st, myData.eventTimes, calcWindow, p.binSize);
trGroupLabels = myData.trGroupLabels;
nGroups = myData.nGroups;
inclRange = bins>p.startRange & bins<=p.stopRange;
spikeCounts = sum(ba(:,inclRange),2)./(p.stopRange-p.startRange);

% PSTH smoothing filter
smWin = p.smWin; 

% smooth ba
baSm = conv2(smWin,1,ba', 'same')'./p.binSize;

% construct psth(s) and smooth it
if p.showAllTraces
    psthSm = zeros(nGroups, numel(bins));
    if p.showErrorShading
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
    if p.showErrorShading
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

if isempty(myData.plotAxes)
    for pidx = 1:3
        subplot(3,1,pidx);
        myData.plotAxes(pidx) = gca;
    end
    set(f, 'UserData', myData);
end

colors = myData.params.colors;
% subplot(3,1,1); 
axes(myData.plotAxes(1));
hold off;
if p.showAllTraces
    for g = 1:nGroups
        plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
        hold on;
    end
else
    plot(bins, psthSm);
end
xlim(plotWindow);
title(['Cluster ' num2str(myData.clusterIDs(p.clusterIndex))]);
xlabel('Time (s)');
ylabel('Firing rate (sp/s)');
yl = ylim();
hold on;
plot(p.startRange*[1 1], yl, 'k--');
plot(p.stopRange*[1 1], yl, 'k--');
makepretty;
box off;

% subplot(3,1,2);
axes(myData.plotAxes(2));
hold off;
plot(rasterX,rasterY, 'k');
xlim(myData.params.window);
ylim([0 length(myData.eventTimes)+1]);
ylabel('event number');
xlabel('time (sec)');
makepretty;
box off;

% subplot(3,1,3);
axes(myData.plotAxes(3));
hold off;
errorbar(trGroupLabels, tuningCurve(:,1), tuningCurve(:,2), 'o-');
xlabel('grouping variable value');
ylabel('average firing rate (Hz)');
makepretty;
box off;

% drawnow;

end

function psthViewerCallback(f, keydata)

% fprintf('callback on %d with source %d\n', f.Number, keydata.Source.Number);



myData = get(f, 'UserData');
p = myData.params;

switch keydata.Key
    case 'rightarrow' % increment cluster index
        
        p.clusterIndex = p.clusterIndex+1;
        if p.clusterIndex>length(myData.clusterIDs)
            p.clusterIndex=1;
        end
        
    case 'leftarrow' % decrement cluster index
        
        p.clusterIndex = p.clusterIndex-1;
        if p.clusterIndex<1
            p.clusterIndex=length(myData.clusterIDs);
        end
        
    case 'uparrow' % increase smoothing
        p.smoothSize = p.smoothSize*1.2;
        p.smWin = genSmWin(p);
        
    case 'downarrow' % decrease smoothing
        p.smoothSize = p.smoothSize/1.2;
        p.smWin = genSmWin(p);
    
    case 'f'        
        p.filterType = p.filterType+1;
        if p.filterType>3; p.filterType = 1; end
        
        p.smWin = genSmWin(p);
        
    case 'e' % whether to show standard error as shading
        p.showErrorShading = ~p.showErrorShading;
        
    case 't' % whether to plot the psth trace for each condition or just the overall one
        p.showAllTraces = ~p.showAllTraces;
        
                    
        
    case 'r'
        ax = subplot(3,1,1); title('click start and stop of range')
        %         [startRange,~] = ginput();
        %         [stopRange,~] = ginput();
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        p.startRange = q(1,1);
        waitforbuttonpress;
        q = get(ax, 'CurrentPoint');
        p.stopRange = q(1,1);
        if p.stopRange<p.startRange
            tmp = p.startRange;
            p.startRange = p.stopRange;
            p.stopRange = tmp;
        end
        
    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            p.clusterIndex = ind;
        end
        
        
end

myData.params = p;
set(f, 'UserData', myData);

% plot with new settings
psthViewerPlot(f)


end



function smWin = genSmWin(p)

switch p.filterType
    case 1 % gaussian
        gw = gausswin(round(p.smoothSize*6),3);
        smWin = gw./sum(gw);
        fprintf(1, 'filter is gaussian with stdev %.2f ms\n', p.smoothSize);
    case 2 % half gaussian, causal
        gw = gausswin(round(p.smoothSize*6),3);
        gw(1:round(numel(gw)/2)) = 0;
        smWin = gw./sum(gw);
        fprintf(1, 'filter is causal half-gaussian with stdev %.2f ms\n', p.smoothSize);
    case 3 % box
        smWin = ones(round(p.smoothSize*2),1);
        fprintf(1, 'filter is box with width %.2f ms\n', p.smoothSize*2);
end
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