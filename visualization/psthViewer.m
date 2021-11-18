

function psthViewer(spikeTimes, clu, eventTimes, window, trGroups, varargin)
% function psthViewer(spikeTimes, clu, eventTimes, window, trGroups[, params])
%
% Controls: see the command window text on startup. 
% params can include 'groupingName' (a string) and 'groupingLegend' (a cell
% array of strings)
%
% TODO:
% XX if plotting all traces, color raster and sort
% - don't replot traces just change X,Y data for speed and for keeping zoom
% levels
% XX indicate on tuning curve which color is which
% XX add ability to switch from graded to distinguishable color scheme, also
% cyclical
% - add support for a second grouping variable - in that case, should use a
% different graded color scheme for each, overlay the traces in the tuning
% curve view, and also provide a 2-d image of tuning curve
% XX add support for plot labels (traces in psth, x-axis in tuning curve)
% XX option to change filter type
% XX expand calc window to avoid falling off at the edges
% XX option for error bars
% XX fix bug where the raster from the last unit stays when there are no new
% spikes to plot

fprintf(1, 'Controls:\n')
fprintf(1, '- left/right arrow: select previous/next cluster\n')
fprintf(1, '- up/down arrow: change smoothing of psth curves\n')
fprintf(1, '- c: dialog box to pick a new cluster ID number\n')
fprintf(1, '- t: toggle showing psth traces for each grouping variable or just the\n')
fprintf(1, 'overall. If showing just overall, raster is sorted chronologically. If\n')
fprintf(1, 'showing by grouping variable, raster is sorted by that variable.\n')
fprintf(1, '- r: select a new range within which to count spikes for the tuning curve\n')
fprintf(1, '- f: change smoothing filter type\n')
fprintf(1, '- x: change color scheme\n')
fprintf(1, '- e: turn on error bars\n')
fprintf(1, '- 1/2: decrease/increase the raster tick size\n')

if ~isempty(varargin)
    p = varargin{1}; % parameters supplied by user
else
    p = [];
end

params.groupingName = getOr(p, 'groupingName', 'Grouping variable value'); 
params.groupingLegend = getOr(p, 'groupingLegend', unique(trGroups)); 

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

params.colorType = 1; 
params.colors = genColors(params.colorType, myData.nGroups); 

myData.params = params;

f = figure; f.Color = 'w';
% set(f, 'Renderer', 'painters')

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
        theseTr = myData.trGroups==trGroupLabels(g);
        psthSm(g,:) = mean(baSm(theseTr,:));
        if p.showErrorShading
            stderr(g,:) = std(baSm(theseTr,:))./sqrt(sum(theseTr));
        end
    end
else
    
    psthSm = mean(baSm);
    if p.showErrorShading
        stderr = std(baSm)./sqrt(size(baSm,1));
    end
    
end

% compute raster
if p.showAllTraces
    [sortedGroups, inds] = sort(myData.trGroups);
    [tr,b] = find(ba(inds,:));
else
    [tr,b] = find(ba);
end
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything

% scale the raster ticks
rasterY(2:3:end) = rasterY(2:3:end)+(myData.params.rasterScale-1);

if isempty(rasterX); rasterX = NaN; rasterY = NaN; end % so that there's something to plot and the plot will clear

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
        if p.showErrorShading
            plotWithErr(bins, psthSm(g,:), stderr(g,:), colors(g,:));
        else
            plot(bins, psthSm(g,:), 'Color', colors(g,:), 'LineWidth', 2.0);
        end
        hold on;
    end
else
    if p.showErrorShading
        plotWithErr(bins, psthSm, stderr, 'k'); 
    else
        plot(bins, psthSm, 'k');
    end
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
if p.showAllTraces
    ug = unique(myData.trGroups);
    
    for g = 1:nGroups
        theseTr = find(sortedGroups==ug(g)); 
        ry = rasterY>=theseTr(1) & rasterY<=theseTr(end);
        
        if any(ry)
            ry(2:3:end) = ry(1:3:end);
            plot(rasterX(ry), rasterY(ry), 'Color', colors(g,:));        
            
        else
            plot(0, 0, 'Color', colors(g,:)); 
        end
        hold on;
    end
else
    plot(rasterX,rasterY, 'k');
end
xlim(myData.params.window);
ylim([0 length(myData.eventTimes)+1]);
ylabel('Event Number');
xlabel('Time (s)');
makepretty;
box off;

% subplot(3,1,3);
axes(myData.plotAxes(3));
hold off;
%errorbar(trGroupLabels, tuningCurve(:,1), tuningCurve(:,2), 'o-');
hh = plot(trGroupLabels, tuningCurve(:,1), 'k'); hold on;
arrayfun(@(x)plot(trGroupLabels(x)*[1 1], tuningCurve(x,1)+tuningCurve(x,2)*[-1 1], 'Color',  colors(x,:), 'LineWidth', 2.0), 1:nGroups);
arrayfun(@(x)plot(trGroupLabels(x), tuningCurve(x,1), 'o', 'Color',  colors(x,:), 'LineWidth', 2.0), 1:nGroups);
xlabel(p.groupingName);
ylabel('Average firing rate (sp/s)');
set(myData.plotAxes(3), 'XTick', trGroupLabels, 'XTickLabel', p.groupingLegend); 
makepretty;
set(hh, 'LineWidth', 1.0); 
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
        
    case 'x'
        p.colorType = p.colorType+1;
        p.colors = genColors(p.colorType, myData.nGroups);
        
    case '2' % increase tick size
        p.rasterScale = p.rasterScale*1.5;
        
    case '1' % decrease tick size
        p.rasterScale = p.rasterScale/1.5;
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
        gw = gausswin(round(p.smoothSize*6*3),3);
        gw(1:round(numel(gw)/2)) = 0;
        smWin = gw./sum(gw);
        fprintf(1, 'filter is causal half-gaussian with stdev %.2f ms\n', p.smoothSize*3);
    case 3 % box
        smWin = ones(round(p.smoothSize*3),1);
        fprintf(1, 'filter is box with width %.2f ms\n', p.smoothSize*3);
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

function colors = genColors(colorType, nColors)

switch mod(colorType, 3)
    case 1 % blue linear map
        colors = copper(nColors); 
        colors = colors(:, [3 2 1]);
    case 2 % distinguishable colors - default matlab order with black and gray
        colors = get(gca, 'ColorOrder'); 
        colors = [colors; 0 0 0; 0.5 0.5 0.5]; 
        nc = size(colors,1);
        colors = colors(mod(0:nColors-1,nc)+1,:);
    case 0 % cyclical 
        %colors = hsv(nColors);
        m = colorcet('C6');
        colors = zeros(nColors,3); 
        for c = 1:3
            qidx = linspace(1,size(m,1), nColors+1); 
            colors(:,c) = interp1(1:size(m,1), m(:,c), qidx(1:end-1));
        end
end
end