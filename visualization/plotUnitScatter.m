
function f = plotUnitScatter(sp, params)
% function plotUnitScatter(sp, params)
% plot of each unit by its basic parameters (amplitude, extent,
% position, firing rate, waveform duration
%
% input argument sp is a struct that's got the kilosort files in it, result
% of loadKSdir
%
% this function needs to be improved - calculations should be taken
% elsewhere and this should just do plotting...
% - also need to double-check that all the indexing is happening correctly
% so things match up 

Fs = sp.sample_rate;

useGood = true;
FRthresh = 0.1;
ampThresh = 75;

if ~isempty(params) && isfield(params, 'useGood')
    useGood = params.useGood;
end

if ~isempty(params) && isfield(params, 'FRthresh')
    FRthresh = params.FRthresh;
end

if ~isempty(params) && isfield(params, 'ampThresh')
    ampThresh = params.ampThresh;
end

% first compute some things
cgs = sp.cgs;
cids = sp.cids;
ycoords = sp.ycoords;

if useGood
    goodCIDs = cids(cgs==2);
    if isempty(goodCIDs)
        % can't use good, there aren't any. Will proceed with all
        goodCIDs = cids;
    end
else
    goodCIDs = cids;
end

clu = sp.clu(ismember(sp.clu, goodCIDs));
spikeTemplates = sp.spikeTemplates(ismember(sp.clu, goodCIDs));

if ~isfield(sp, 'tempPerClu')
    tempPerClu = findTempForEachClu(clu, spikeTemplates);    
else
    tempPerClu = sp.tempPerClu;
end

if ~isfield(sp, 'cluYpos')
    sp.cluYpos = sp.templateYpos(tempPerClu+1);
end
if ~isfield(sp, 'cluAmps')
    sp.cluAmps = sp.tempAmps(tempPerClu+1);
end

temps = sp.tempsUnW;
tempChanAmps = squeeze(max(temps,[],2))-squeeze(min(temps,[],2));

nChansSpan = zeros(1, length(goodCIDs));
for c = 1:length(goodCIDs)
    thisTempInd = tempPerClu(goodCIDs(c)+1);
    thisTempAmps = tempChanAmps(thisTempInd+1,:);
    ampMed = median(thisTempAmps);
    peakYC = ycoords(thisTempAmps==max(thisTempAmps));
    peakYC = peakYC(1);
    thisTempAmps(ycoords<peakYC-150 | ycoords>peakYC+150) = 0;
    
    nChansSpan(c) = sum(thisTempAmps>6*ampMed);
end
sp.nChansSpan = nChansSpan;

if ~isfield(sp, 'FRs')
    [cids, spikeCounts] = countUnique(clu);
    sp.FRs = spikeCounts/max(sp.st); % convert to Hz
end

goodAmps = sp.cluAmps(ismember(cids, goodCIDs));
goodFRs = sp.FRs(ismember(cids, goodCIDs));
goodYpos = sp.cluYpos(ismember(cids, goodCIDs));
nChansSpan = sp.nChansSpan;
goodDur = sp.tempDur(tempPerClu(goodCIDs+1)+1);
isNarrow = goodDur/Fs<0.00033; % a third of a ms

inclUnits = goodAmps>ampThresh & goodFRs>FRthresh;

goodAmps = goodAmps(inclUnits);
goodFRs = goodFRs(inclUnits);
goodYpos = goodYpos(inclUnits);
nChansSpan = nChansSpan(inclUnits);
isNarrow = isNarrow(inclUnits);

f = figure;

spanScaleFactor = 4;

% these plots are for the legend
h = scatter(-1000, -1000, spanScaleFactor, 50, 'filled'); hold on;
set(h, 'MarkerEdgeColor', 'k');
h = scatter(-1100, -1100, spanScaleFactor*5, 50, 'filled'); hold on;
set(h, 'MarkerEdgeColor', 'k');
h = scatter(-1200, -1200, spanScaleFactor*15, 50, 'filled'); hold on;
set(h, 'MarkerEdgeColor', 'k');
h = scatter(-1300, -1300, spanScaleFactor*5, 50, '^', 'filled'); hold on;
set(h, 'MarkerEdgeColor', 'k');
h = scatter(-1400, -1400, spanScaleFactor*5, 50, 'filled'); hold on;
set(h, 'MarkerEdgeColor', 'k');

h = scatter(goodAmps(~isNarrow), goodYpos(~isNarrow), nChansSpan(~isNarrow)*spanScaleFactor, goodFRs(~isNarrow), '^', 'filled');
set(h, 'MarkerEdgeColor', 'k');

hold on
h = scatter(goodAmps(isNarrow), goodYpos(isNarrow), nChansSpan(isNarrow)*spanScaleFactor, goodFRs(isNarrow), 'filled');
set(h, 'MarkerEdgeColor', 'k');

ylim([0 3840])
caxis([0 40]);
colormap hot
h = colorbar;
h.Label.String = 'firing rate (sp/s)';
xlabel('waveform amplitude (µV)');
[h,icons,plots,legend_text] = legend({'small span', 'medium span', 'large span', 'broad waveform', 'narrow waveform'}, 'Location', 'NorthWestOutside');
icons(6).Children.MarkerSize = icons(6).Children.MarkerSize/2;
icons(8).Children.MarkerSize = icons(8).Children.MarkerSize*2;
makepretty;
