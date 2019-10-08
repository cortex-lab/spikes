

function f = evRastersGUI(st, clu, cweA, cwtA, moveData, lickTimes, anatData)
% function evRastersGUI(st, clu, cweA, cwtA, moveData, lickTimes, anatData)
%
% Displays rasters and PSTHs for individual units in the choiceworld task
%
% Inputs:
% - st - vector of spike times
% - clu - vector of cluster identities
% - cweA - table of trial labels, containing contrastLeft, contrastRight,
% choice, and feedback
% - cwtA - table of times of events in trials, containing stimOn, beeps,
% and feedbackTime
% - moveData - a struct with moveOnsets, moveOffsets, moveType
% - lickTimes - a vector of lick times
% - anatData - a struct with: 
%   - coords - [nCh 2] coordinates of sites on the probe
%   - wfLoc - [nClu nCh] size of the neuron on each channel
%   - borders - table containing upperBorder, lowerBorder, acronym
%   - clusterIDs - an ordering of clusterIDs that you like (wfLoc and
%   waveforms must be sorted to match this order!!)
%   - waveforms - [nClu nCh nTimepoints] waveforms of the neurons
%
% Controls: 
% - up/down arrows to switch between clusters
% - left/right arrows to increase/decrease smoothing size
% - 'h' to hide/unhide behavioral icons on the rasters
% - t/r to increase/decrease raster tick sizes
% - 'c' to go to a particular cluster by number
% - click on anatomy to jump to cluster nearest there
% - 'l' to see the legend

% TODO
% - click on stim onset plots to set counting window for contrast resp
% function

% fprintf(1, 'Controls: \n');
% fprintf(1, ' - up/down arrows to switch between clusters\n');
% fprintf(1, ' - left/right arrows to increase/decrease smoothing size\n');
% fprintf(1, ' - ''h'' to hide/unhide behavioral icons on the rasters\n');
% fprintf(1, ' - t/r to increase/decrease raster tick sizes\n');
% fprintf(1, ' - ''c'' to go to a particular cluster by number\n');
% fprintf(1, ' - click on anatomy to jump to cluster nearest there\n');


% construct figure
f = figure; set(f, 'Color', 'w');


pars.nTop = 6; % num subplots in the top two rows
pars.nBottom = 7;
pars.nVertSp = 6; % the four rows get 2, 1, 2, 1
pars.psthBinSize = 0.002;
pars.smoothWinStd = 0.03;
% pars.smoothWin = myGaussWin(pars.smoothWinStd, 1/pars.psthBinSize);
pars.smoothWin = createSmWin(pars.smoothWinStd, pars.psthBinSize);
pars.lickBoutGap = 0.25;
pars.evsVisible = true;
pars.cluIndex = 1;
pars.nWFtoPlot = 24;
pars.crfWin = [0.025 0.225]; % sec, window around stimOnset to compute contrast response function

myData.st = st;
myData.clu = clu;
if ~isempty(anatData) && isfield(anatData, 'clusterIDs')
    myData.clusterIDs = anatData.clusterIDs;
else
    myData.clusterIDs = unique(clu);
end
myData.cweA = cweA;
myData.cwtA = cwtA;
pars.tickSize = ceil(numel(cwtA.stimOn)/300);
myData.moveOnsets = moveData.moveOnsets(:);
myData.moveOffsets = moveData.moveOffsets(:);
moveData.moveType(moveData.moveType==3) = 0; % call the unclassified ones flinches
myData.moveType = moveData.moveType(:);
myData.lickTimes = lickTimes;
myData.anatData = anatData;
myData.pars = pars;
myData.f = f;

ev = createEv(myData);
myData.ev = ev;

evData = createEvTimesOrders(myData);
myData.evData = evData; 

myData = createPlots(myData);

% precompute all the BAs here?

updatePlots(myData);

set(f, 'UserData', myData);
set(f, 'KeyPressFcn', @(f,k)kpCallback(f, k));

function kpCallback(f,keydata)
myData = get(f, 'UserData');
switch keydata.Key
    case 'uparrow' % increment cluster index
        
        myData.pars.cluIndex = myData.pars.cluIndex+1;
        if myData.pars.cluIndex>length(myData.clusterIDs)
            myData.pars.cluIndex=1;
        end        
        updatePlots(myData);

    case 'downarrow' % decrement cluster index
        
        myData.pars.cluIndex = myData.pars.cluIndex-1;
        if myData.pars.cluIndex<1
            myData.pars.cluIndex=length(myData.clusterIDs);
        end
        updatePlots(myData);

    case 'rightarrow' % wider smoothing
        
        myData.pars.smoothWinStd = myData.pars.smoothWinStd*5/4;
        %myData.pars.smoothWin = myGaussWin(myData.pars.smoothWinStd, 1/myData.pars.psthBinSize);     
        myData.pars.smoothWin = createSmWin(myData.pars.smoothWinStd, myData.pars.psthBinSize);
        updatePlots(myData);

    case 'leftarrow' % narrower smoothing
        
        myData.pars.smoothWinStd = myData.pars.smoothWinStd*4/5;
        %myData.pars.smoothWin = myGaussWin(myData.pars.smoothWinStd, 1/myData.pars.psthBinSize);     
        myData.pars.smoothWin = createSmWin(myData.pars.smoothWinStd, myData.pars.psthBinSize);
        updatePlots(myData);

    case 'h'
        myData.pars.evsVisible = ~myData.pars.evsVisible;
        for p = 1:numel(myData.hRasterEvs)
            for q = 1:numel(myData.hRasterEvs{p})
                if myData.pars.evsVisible
                    set(myData.hRasterEvs{p}(q), 'Visible', 'on');
                else
                    set(myData.hRasterEvs{p}(q), 'Visible', 'off');
                end
            end
        end
        
    case 't' % bigger raster ticks
        myData.pars.tickSize = myData.pars.tickSize*2;
        updatePlots(myData);

    case 'r' % smaller raster ticks
        myData.pars.tickSize = myData.pars.tickSize/2;
        updatePlots(myData);
        
    case 'c'
        newC = inputdlg('cluster ID?');
        ind = find(myData.clusterIDs==str2num(newC{1}),1);
        if ~isempty(ind)
            myData.pars.cluIndex = ind;
        else
            fprintf(1, 'could not find cluster %d\n', str2num(newC{1}));
        end
        updatePlots(myData);

    case 'l'
        ev = myData.ev;
        f = figure; f.Name = 'Behavioral Icons legend';
        f.ToolBar = 'none'; f.MenuBar = 'none';
        for e = 1:numel(ev)
            plot(0,0,ev(e).icon, 'Color', ev(e).color); 
            hold on; 
        end
        legend({ev.name});
        axis off;
        
    case 'u'
        updatePlots(myData);
end

% updatePlots(myData);
set(f, 'Name', sprintf('clusterID = %d', myData.clusterIDs(myData.pars.cluIndex)));

set(f, 'UserData', myData);

function myData = createPlots(myData)
nTop = myData.pars.nTop; % num subplots in the top two rows
nBottom = myData.pars.nBottom;
nVertSp = myData.pars.nVertSp;

for colInd = 1:nTop
    axRaster(colInd) = subplot(nVertSp, nTop, [colInd nTop+colInd]);
    axPSTH(colInd) = subplot(nVertSp, nTop, nTop*2+colInd);
end
for colInd = 1:nBottom
    axRaster(colInd+nTop) = subplot(nVertSp, nBottom, nBottom*[3 4]+colInd);
    axPSTH(colInd+nTop) = subplot(nVertSp, nBottom, nBottom*5+colInd);
end

% shift them all to the right just a bit to make space for waveform on the
% left
for n = 1:numel(axRaster)
    p = get(axRaster(n), 'Position'); 
    set(axRaster(n), 'Position', [p(1)+0.07 p(2) p(3) p(4)]);
    p = get(axPSTH(n), 'Position'); 
    set(axPSTH(n), 'Position', [p(1)+0.07 p(2) p(3) p(4)]);
end

evData = myData.evData;
for e = 1:length(evData)
    hRaster(e) = plot(axRaster(e), 0, 0, 'k');
    hold(axRaster(e), 'on');
    hRasterEvs{e} = addEvents(axRaster(e), evData(e).times, evData(e).trOrders, myData.ev, evData(e).windows);       
    plot(axRaster(e), [0 0], [0 numel(evData(e).times)], 'k');
    ylim(axRaster(e), [0 numel(evData(e).times)]+0.5);
    xlim(axRaster(e), evData(e).windows);
    box(axRaster(e), 'off');   
    ylabel(axRaster(e), evData(e).ylab);
end
myData.hRaster = hRaster;
myData.hRasterEvs = hRasterEvs;

%  create psth traces
for e = 1:length(evData)
    hold(axPSTH(e), 'on');
    for c = 1:size(evData(e).colors)
        hPSTH{e}(c) = plot(axPSTH(e), evData(e).windows, [0 0], ...
            'Color', evData(e).colors(c,:), 'LineWidth',2.0);    
        if ~isempty(evData(e).lineType)
            hPSTH{e}(c).LineStyle = evData(e).lineType{c};
        end
    end
    xlim(axPSTH(e), evData(e).windows);
    box(axPSTH(e), 'off');
    hold(axPSTH(e), 'on');
    plot(axPSTH(e), [0 0], [0 1000], 'k'); %1000 is supposed to just be bigger than it ever could - axes will rescale from unit to unit
    xlabel(axPSTH(e), sprintf('time from %s (sec)', evData(e).alignName));
    
end
myData.hPSTH = hPSTH;
myData.axPSTH = axPSTH;

% contrast response function
axCRF = axes(myData.f);
set(axCRF, 'Position', [0.28 0.2 0.05 0.1]);
imCRF = imagesc(axCRF, rand(4));
colormap(axCRF, parula);
axis(axCRF, 'off');
set(axCRF, 'YDir', 'normal');
myData.axCRF = axCRF;
myData.imCRF = imCRF;

if isfield(myData, 'anatData')
    anatData = myData.anatData;
    axAnat = axes(myData.f);
    set(axAnat, 'Position', [0.01 0.05 0.03 0.9]);
    
    wfSize = 0.1*ones(size(anatData.coords(:,1))); 
    hProbeScatter = scatter(anatData.coords(:,1), anatData.coords(:,2), wfSize);
    hProbeScatter.MarkerFaceColor = 'flat';
    set(hProbeScatter, 'HitTest', 'off');
    
    hold on; 
    for b = 1:size(anatData.borders, 1)
        h = plot([min(anatData.coords(:,1)) max(anatData.coords(:,1))], ...
            anatData.borders.upperBorder(b)*[1 1], 'k');
        set(h, 'HitTest', 'off');
        h = plot([min(anatData.coords(:,1)) max(anatData.coords(:,1))], ...
            anatData.borders.lowerBorder(b)*[1 1], 'k');
        set(h, 'HitTest', 'off');
        ah = annotation('textbox', 'String', anatData.borders.acronym{b});
        ah.EdgeColor = 'none';
        set(ah, 'Parent', axAnat);
        set(ah, 'Position', [max(anatData.coords(:,1))+10, mean([anatData.borders.lowerBorder(b), anatData.borders.upperBorder(b)]), 0.05, 0.05])
        set(ah, 'HitTest', 'off');
    end    
    %axis(axAnat, 'off')
    drawnow;
    axAnat.XRuler.Axle.LineStyle = 'none'; 
    axAnat.YRuler.Axle.LineStyle = 'none'; 
    set(axAnat, 'YTick', [], 'XTick', []);
    set(axAnat, 'ButtonDownFcn', @(q,k)anatClick(q, k, myData.f));
    
    myData.hProbeScatter = hProbeScatter;
    myData.axAnat = axAnat;
    
    if isfield(anatData, 'waveforms')
        axWF = axes(myData.f);
        set(axWF, 'Position', [0.07 0.55 0.1 0.4]);
        hold(axWF, 'on');
        for n = 1:myData.pars.nWFtoPlot
            hWF(n) = plot(axWF, 0, 0, 'k', 'LineWidth', 2.0);
        end
        axis(axWF, 'off');
        myData.hWF = hWF;
        
        
        axWFslice = axes(myData.f);
        set(axWFslice, 'Position', [0.07 0.40 0.1 0.15]);
        hSliceIm = imagesc(axWFslice, rand(10));
        colormap(axWFslice, colormap_BlueWhiteRed); axis(axWFslice, 'off');
        hold(axWFslice, 'on');
        set(axWFslice, 'YDir','normal');
        plot(axWFslice, size(anatData.waveforms,3)*[1 1], [0 1000], 'k');
        myData.axWFslice = axWFslice;
        myData.hSliceIm = hSliceIm;
        
    end
end
axACGlin = axes(myData.f);
set(axACGlin, 'Position', [0.07 0.23 0.1 0.15]);
axACGlog = axes(myData.f);
set(axACGlog, 'Position', [0.07 0.04 0.1 0.15]);
myACG(rand(1,100), axACGlin, axACGlog);
myData.axACGlin = axACGlin;
myData.axACGlog = axACGlog;



function anatClick(~, keydata, f)

clickY = keydata.IntersectionPoint(2);
myData = get(f, 'UserData');

% find the cluster with the closest position to that click location
wfLoc = myData.anatData.wfLoc;
coords = myData.anatData.coords;
[~, maxChan] = max(wfLoc, [], 2);
maxY = coords(maxChan,2);
[~, closestInd] = min(abs(clickY-maxY));
myData.pars.cluIndex = closestInd;
updatePlots(myData);
set(f, 'Name', sprintf('clusterID = %d', myData.clusterIDs(myData.pars.cluIndex)));
set(f, 'UserData', myData);

function updatePlots(myData)
% pick spike times and compute the BAs
st = myData.st(myData.clu==myData.clusterIDs(myData.pars.cluIndex));

[allBA, allBins] = computeBAs(st, myData);


evData = myData.evData;
maxyl = [1000 -1000];
for e = 1:length(evData)
    ba = allBA{e}; bins = allBins{e};
    ba = ba(evData(e).trOrders,:);
    
    % set the new rasters
    [tr,b] = find(ba);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy*myData.pars.tickSize+reshape(repmat(tr',3,1),1,length(tr)*3)-0.5; % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
    set(myData.hRaster(e), 'XData', rasterX, 'YData', rasterY);
    
    % set the new PSTH traces
    %gIDs = unique(evData(e).groups); gIDs = gIDs(~isnan(gIDs));
    gIDs = evData(e).gIDs;
    for g = 1:length(gIDs)
        inclTrials = ~isnan(evData(e).times)&evData(e).groups==gIDs(g);
        if ~isempty(inclTrials)
            thisPSTH = conv(nanmean(ba(inclTrials,:),1)./myData.pars.psthBinSize, myData.pars.smoothWin, 'same');
            set(myData.hPSTH{e}(g), 'XData', bins, 'YData', thisPSTH);
            if max(thisPSTH)>maxyl(2)
                maxyl(2) = max(thisPSTH);
            end
            if min(thisPSTH)<maxyl(1)
                maxyl(1) = min(thisPSTH);
            end
        end
    end    
end
if maxyl(2)>maxyl(1) % can fail when there were no spikes in any window
    for e = 1:length(evData)
        ylim(myData.axPSTH(e), maxyl);
    end
end

% now the contrast response function. Here we're not supposed to know
% anything about the particulars of the data in ba and evData, but this is
% a special case I guess. 

ba = allBA{2}; bins = allBins{2}; % stim-aligned
cR = myData.cweA.contrastRight; cL = myData.cweA.contrastLeft;
uR = unique(cR); nR = numel(uR);
uL = unique(cL); nL = numel(uL);
crf = nan(nR, nL); 
inclBins = bins>myData.pars.crfWin(1) & bins<myData.pars.crfWin(2);
for ri = 1:nR
    for li = 1:nL
        crf(li,ri) = mean(mean(ba(cR==uR(ri)&cL==uL(li),inclBins)));
    end
end
set(myData.imCRF, 'CData', crf);

if isfield(myData, 'hProbeScatter')
    anatData = myData.anatData;
    minSize = 0.1; maxSize = 30;
    wfSize = anatData.wfLoc(myData.pars.cluIndex,:);
    % normalize to range
    wfSize = (wfSize-min(wfSize))/(max(wfSize)-min(wfSize))*(maxSize-minSize)+minSize;
    myData.hProbeScatter.SizeData = wfSize;
end

% new acg
[xLin, nLin, ~, nLog] = myACG(st,[],[]);
ch = get(myData.axACGlin, 'Children');
nLin = nLin(xLin<0.1);
yyLin = reshape([nLin nLin]', numel(nLin)*2,1);
set(ch(3), 'YData', yyLin(1:end-1));
set(ch(1), 'YData', nLog(end)*[1 1]);
set(ch(2), 'YData', [0 max(yyLin)]);
if max(yyLin)>0; set(myData.axACGlin, 'YLim', [0 max(yyLin)]); end

ch = get(myData.axACGlog, 'Children');
yyLog = reshape([nLog nLog]', numel(nLog)*2,1);
set(ch(3), 'YData', yyLog);
set(ch(1), 'YData', nLog(end)*[1 1]);
set(ch(2), 'YData', [0 max(yyLog)]);
if max(yyLog)>0; set(myData.axACGlog, 'YLim', [0 max(yyLog)]); end

if isfield(myData, 'hWF')
    % determine nearest channels
    wfSize = anatData.wfLoc(myData.pars.cluIndex,:);
    [~,maxCh] = max(wfSize);
    c = anatData.coords;
    chDists = ((c(maxCh,1)-c(:,1)).^2+(c(maxCh,2)-c(:,2)).^2).^0.5;
    [~,nearestCh] = sort(chDists);
    % update the traces with these wfs
    thisWF = squeeze(anatData.waveforms(myData.pars.cluIndex,:,:));
    for n = 1:myData.pars.nWFtoPlot
        set(myData.hWF(n), 'YData', thisWF(nearestCh(n),:)*5+c(nearestCh(n),2), ...
            'XData', (0:size(thisWF,2)-1)*0.45+c(nearestCh(n),1));
    end
    
    % new wf slice
    wfSlice1 = waveformSpaceTimeSlice(thisWF, c, 19);
    wfSlice2 = waveformSpaceTimeSlice(thisWF, c, 51);
    wfSlice = [wfSlice1 wfSlice2];
    set(myData.hSliceIm, 'CData', wfSlice);
    xlim(myData.axWFslice,[0 size(wfSlice,2)+1])
    ylim(myData.axWFslice,[0 size(wfSlice,1)+1])
    caxis(myData.axWFslice, max(abs(wfSlice(:)))*0.2*[-1 1]);
end
function [ba, bins] = computeBAs(st, myData)
% compute all the binned arrays of spikes for new spike times
evData = myData.evData;
for e = 1:length(evData)
    eventTimes = evData(e).times;
    
    % is this one that we already computed? 
    % Note this has a problem if you want to use a bigger window the second
    % time around. For now, ignoring this. 
    found = false;
    for e2 = 1:e-1
        if ~found && numel(eventTimes)==numel(evData(e2).times) && ...
                sum(eventTimes==evData(e2).times)==numel(eventTimes)
            ba{e} = ba{e2};
            bins{e} = bins{e2};
            found = true;
        end
    end
       
    if ~found
        nTimes = length(eventTimes);
        windowExp = evData(e).windows+myData.pars.smoothWinStd*5*[-1 1];
        [newBA, newBins] = timestampsToBinned(st, eventTimes, myData.pars.psthBinSize, windowExp);
        ba{e} = newBA;
        bins{e} = newBins;
    end
end

function h = addEvents(ax, eventTimes, trOrder, ev, thisWindow)
% add other events nearby
hold(ax, 'on');
reOrderEvent = eventTimes(trOrder);
nTimes = numel(reOrderEvent);
hInd = 1;
for e = 1:length(ev)
    otherEvent = ev(e).times;
    x = WithinRanges(otherEvent, bsxfun(@plus, reOrderEvent, thisWindow), [1:nTimes]');
    [ii,trInds] = find(x);
    relTimes = otherEvent(ii)-reOrderEvent(trInds);
    q = plot(ax, relTimes, trInds, ev(e).icon, 'Color', ev(e).color);
    if ~isempty(q)
        h(hInd) = q;
        hInd = hInd+1;
    end
end

function evData = createEvTimesOrders(myData)
cweA = myData.cweA; cwtA = myData.cwtA;
moveOnsets = myData.moveOnsets(:); moveType = myData.moveType(:); moveOffsets = myData.moveOffsets(:);
lickTimes = myData.lickTimes;

% trial ordering for top row
stimOn = cwtA.stimOn;
respMoves = moveOnsets(moveType==1 | moveType==2);
nextMoveTime = arrayfun(@(x)respMoves(find(respMoves>x,1)), stimOn);
[~,trOrder] = sort(nextMoveTime-stimOn);

n = 0;

% top row, move offset before trial start
n = n+1;
if moveOffsets(1)>stimOn(1); moveOffsets = [0; moveOffsets]; end
prevMoveTime = arrayfun(@(x)moveOffsets(find(moveOffsets<x,1, 'last')), stimOn);
evData(n).times = prevMoveTime;
evData(n).trOrders = [1:numel(prevMoveTime)]';
evData(n).windows = [-0.3 0.7];
evData(n).groups = ones(size(trOrder));
evData(n).gIDs = 1;
evData(n).colors = [0 0 0];
evData(n).alignName = 'pre-trial move offset';
evData(n).ylab = 'trial chronologically';


% top row, stimulus
n = n+1;
evData(n).times = cwtA.stimOn;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = ones(size(trOrder));
evData(n).gIDs = 1;
evData(n).colors = [0 0 0];
evData(n).alignName = 'stim onset';
evData(n).ylab = 'trial by (moveOnset-stimOn)';

% top row, first movement
stimOn = cwtA.stimOn;
endTime = cwtA.feedbackTime;
nextMoveInd = arrayfun(@(x)find(moveOnsets>x,1), stimOn);
nextMoveTime = moveOnsets(nextMoveInd);
nextMoveTime(nextMoveTime>endTime) = NaN;
n = n+1;
evData(n).times = nextMoveTime;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.5 0.5];
evData(n).groups = ones(size(trOrder));
evData(n).gIDs = 1;
evData(n).colors = [0 0 0];
evData(n).alignName = 'first move';
evData(n).ylab = 'trial by (moveOnset-stimOn)';

% top row, go cue
n = n+1;
evData(n).times = cwtA.beeps;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = ones(size(trOrder));
evData(n).gIDs = 1;
evData(n).colors = [0 0 0];
evData(n).alignName = 'go cue';
evData(n).ylab = 'trial by (moveOnset-stimOn)';

% top row, feedback
n = n+1;
evData(n).times = cwtA.feedbackTime;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = ones(size(trOrder));
evData(n).gIDs = 1;
evData(n).colors = [0 0 0];
evData(n).alignName = 'feedback';
evData(n).ylab = 'trial by (moveOnset-stimOn)';

% top row, licks
lickBoutGap = myData.pars.lickBoutGap; % seconds, if gap is <this then consider it part of the same bout
n = n+1;
evData(n).times = myData.lickTimes(diff([0; myData.lickTimes])>lickBoutGap);
evData(n).trOrders = [1:numel(evData(n).times)]';
evData(n).windows = [-0.3 0.7];
evData(n).groups = ones(size(evData(n).trOrders));
evData(n).gIDs = 1;
evData(n).colors = [0 0 0];
evData(n).alignName = 'licks';
evData(n).ylab = 'licks (chronologically)';

% bottom row: stimulus left
[sortedContrasts,trOrder] = sortrows([cweA.contrastLeft cweA.contrastRight], [1 2]);
n = n+1;
groups = nan(size(trOrder));
cL = unique(cweA.contrastLeft); cR = unique(cweA.contrastRight);
for con = 1:4
    inclTrials = ~isnan(stimOn)&sortedContrasts(:,1)==cL(con)&sortedContrasts(:,2)==0;
    groups(inclTrials)=con;
end
evData(n).times = cwtA.stimOn;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = groups;
evData(n).gIDs = [1 2 3 4];
visColorsL = copper(4); visColorsL = visColorsL(2:4, [3 1 2]);
evData(n).colors = [0 0 0; visColorsL];
evData(n).alignName = 'stim onset LEFT';
evData(n).ylab = 'trial by stim left contrast';

% bottom row: stimulus right
[sortedContrasts,trOrder] = sortrows([cweA.contrastLeft cweA.contrastRight], [2 1]);
n = n+1;
groups = nan(size(trOrder));
for con = 1:4
    inclTrials = ~isnan(stimOn)&sortedContrasts(:,2)==cR(con)&sortedContrasts(:,1)==0;
    groups(inclTrials)=con;
end
evData(n).times = cwtA.stimOn;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = groups;
evData(n).gIDs = [1 2 3 4];
visColorsR = copper(4); visColorsR = visColorsR(2:4, [1 3 2]);
evData(n).colors = [0 0 0; visColorsR];
evData(n).alignName = 'stim onset RIGHT';
evData(n).ylab = 'trial by stim right contrast';


% bottom row: equal contrasts, aligned to stim, split by choice
incl = cweA.contrastLeft==cweA.contrastRight;
con = cweA.contrastLeft(incl);
choice = cweA.choice(incl);
[sortedTr,trOrder] = sortrows([choice con]);
n = n+1;
evData(n).times = cwtA.stimOn(incl);
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = sortedTr(:,1);
evData(n).gIDs = [1 2 3];
evData(n).colors = [0 0.5 1; 1 0.5 0; 0 0 0];
evData(n).alignName = 'stim onset EQ CON';
evData(n).ylab = 'eq con trial by choice (L, R, nogo)';

% bottom row: equal contrasts, aligned to movement onset, split by choice
stimOn = cwtA.stimOn;
endTime = cwtA.feedbackTime;
nextMoveInd = arrayfun(@(x)find(moveOnsets>x,1), stimOn);
nextMoveTime = moveOnsets(nextMoveInd);
nextMoveTime(nextMoveTime>endTime) = NaN;
nextMoveType = moveType(nextMoveInd); 

incl = cweA.contrastLeft==cweA.contrastRight & ~isnan(nextMoveTime);
[sortedTr,trOrder] = sort(nextMoveType(incl));
n = n+1;
evData(n).times = nextMoveTime(incl);
evData(n).trOrders = trOrder;
evData(n).windows = [-0.5 0.5];
evData(n).groups = sortedTr;
evData(n).gIDs = [0 1 2];
evData(n).colors = [0 0 0; 0 0.5 1; 1 0.5 0];
evData(n).alignName = 'first move EQ CON';
evData(n).ylab = 'eq con trial by first move type (flinch, L, R)';


% bottom row: move onsets
movesInTrials = logical(WithinRanges(moveOnsets, [cwtA.stimOn cwtA.feedbackTime]));
theseMoveTimes = moveOnsets(movesInTrials&moveType<3);
theseMoveTypes = moveType(movesInTrials&moveType<3);
[sortedTypes,trOrder] = sort(theseMoveTypes);
n = n+1;
evData(n).times = theseMoveTimes;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.5 0.5];
evData(n).groups = sortedTypes;
evData(n).gIDs = [0 1 2];
evData(n).colors = [0 0 0; 0 0.5 1; 1 0.5 0];
evData(n).alignName = 'moves within trial';
evData(n).ylab = 'trial by move type (flinch, L, R)';


% bottom row: go cues
choice = cweA.choice;
[sortedChoices,trOrder] = sort(choice);
n = n+1;
evData(n).times = cwtA.beeps;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = sortedChoices;
evData(n).gIDs = [1 2 3];
evData(n).colors = [0 0.5 1; 1 0.5 0; 0 0 0];
evData(n).alignName = 'go cue';
evData(n).ylab = 'trial by choice (L, R, nogo)';


% bottom row: rewards
feedbackType = cweA.feedback;
[sortedFeedback,trOrder] = sortrows([choice feedbackType]);
n = n+1;
evData(n).times = cwtA.feedbackTime;
evData(n).trOrders = trOrder;
evData(n).windows = [-0.3 0.7];
evData(n).groups = sortedFeedback(:,1).*sortedFeedback(:,2); %this works because feedback is -1 or 1 so these come out uniquely
evData(n).gIDs = [-3 -2 -1 1 2 3];
% evData(n).colors = [1 0.3 0.3; 0 0.8 0.5];
evData(n).colors = [0 0 0; 1 0.5 0; 0 0.5 1; 0 0.5 1; 1 0.5 0; 0 0 0];
evData(n).lineType = {':', ':', ':', '-','-','-'};
evData(n).alignName = 'feedback';
evData(n).ylab = 'trial by feedback type (neg, reward)';




function ev = createEv(myData)
cweA = myData.cweA; cwtA = myData.cwtA;
moveOnsets = myData.moveOnsets; moveType = myData.moveType;
lickTimes = myData.lickTimes;
% Icons and times for events
% visual stimuli: plus for left, circles for right, so you can see when
% there's both
cL = unique(cweA.contrastLeft); cR = unique(cweA.contrastRight);
visColorsL = copper(4); visColorsL = visColorsL(2:4, [3 1 2]);
for con = 1:3
    ev(con).name = sprintf('visL%d', con);
    ev(con).icon = '*'; 
    ev(con).color = visColorsL(con,:);
    ev(con).times = cwtA.stimOn(cweA.contrastLeft==cL(con+1));
end
visColorsR = copper(4); visColorsR = visColorsR(2:4, [1 3 2]);
for con = 1:3
    ev(con+3).name = sprintf('visR%d', con);
    ev(con+3).icon = 'o'; 
    ev(con+3).color = visColorsR(con,:);
    ev(con+3).times = cwtA.stimOn(cweA.contrastRight==cR(con+1));
end

moveTypes = {'flinch', 'moveL', 'moveR', 'other'};
moveColors = [0 0 0; 0 0.5 1; 1 0.5 0; 0.8 0.8 0.8]; 
for m = 1:4
    n = length(ev)+1;
    ev(n).name = moveTypes{m};
    ev(n).icon = '^';
    ev(n).color = moveColors(m,:);
    ev(n).times = moveOnsets(moveType==m-1);
end

goColors = [0 0.5 1; 1 0.5 0; 0 0 0]; 
for r = 1:3
    n = length(ev)+1;
    ev(n).name = sprintf('choice%d', r);
    ev(n).icon = 'x';
    ev(n).color = goColors(r,:);
    ev(n).times = cwtA.beeps(cweA.choice==r);
end

rewardIcon = 'p';
rewardColors = [0 0.8 0.5; 1 0.3 0.3];

n = length(ev)+1;
ev(n).name = 'reward';
ev(n).icon = rewardIcon;
ev(n).color = rewardColors(1,:);
ev(n).times = cwtA.feedbackTime(cweA.feedback==1);
n = length(ev)+1;
ev(n).name = 'negFeedback';
ev(n).icon = rewardIcon;
ev(n).color = rewardColors(2,:);
ev(n).times = cwtA.feedbackTime(cweA.feedback==-1);

lickBoutGap = myData.pars.lickBoutGap; % seconds, if gap is <this then consider it part of the same bout
lickBoutStarts = lickTimes(diff([0; lickTimes])>lickBoutGap);
n = length(ev)+1;
ev(n).name = 'lickBout';
ev(n).icon = '.';
ev(n).color = [1 0 1];
ev(n).times = lickBoutStarts;


function smWin = createSmWin(stdev, psthBinSize)
smWin = myGaussWin(stdev, 1/psthBinSize);
smWin(1:round(numel(smWin)/2)) = 0;
smWin = smWin./sum(smWin);
