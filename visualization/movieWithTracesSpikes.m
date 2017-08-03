
% todo
% - implement pixel traces
% - implement version with raster
% - smoothing
% - normalization 
% - implement 'g' for go to specific time

function movieWithTracesSpikes(spikeTimes, spikeClu, traces, auxVid, anatData, pars)
% function movieWithTracesSpikes(spikeTimes, spikeClu, traces, auxVid, pars)
% - spikeTimes is nSpikes x 1
% - spikeClu is nSpikes x 1
% - traces is a struct array with three fields:
%   - t is 1 x nTimePoints
%   - v is 1 x nTimePoints
%   - name is a string
% - auxVid is a struct array with fields:
%   - data is a cell with any required data, e.g. {readerObject, timestamps}
%   - f is a function handle that will be called to show a frame like this:
%       f(ax, currentTime, data)
%   - name is a string
%
% Usage:
% - 'p' plays/pauses
% - 'r' starts/stops recording, if playing
% - up/down arrow keys increase/decrease playback rate
% - alt+arrowkeys changes view
% - click to change the location of the plotted point. 
%   - right click to add a new point
%   - 'c' to clear the plotted points, leaving only the last one
% - '-' and '=' change the caxis, scaling up and down. It will stay
% centered around zero though. 
% - 'b' jumps backwards

if isempty(pars); pars = struct(); end;

ud.binSpikes = pick(pars, 'binSpikes', 'def', true);
ud.binSize = pick(pars, 'binSize', 'def', 0.02);
movieSaveFilePath = pick(pars, 'movieSaveFilePath', 'def', []);
ud.winSize = pick(pars, 'winSize', 'def', [-5 1]);
ud.nBinInWin = round(diff(ud.winSize)/ud.binSize);
startTime = pick(pars, 'startTime', 'def', ud.binSize-ud.winSize(1));
ud.normalize = pick(pars, 'normalize', 'def', true);
ud.smoothSizeT = pick(pars, 'smoothSizeT', 'def', 0);

allData.st = spikeTimes;
allData.clu = spikeClu;

if ud.binSpikes
    fprintf(1, 'binning spikes...\n');
    allData.binnedSpikes = full(sparse(ceil(allData.st./ud.binSize), spikeClu, ones(size(allData.st))));
    allData.tInds = [0:size(allData.binnedSpikes,1)-1]*ud.binSize;
    
    if ud.normalize
        fprintf(1, 'normalizing...\n');
        mn = mean(allData.binnedSpikes);
%         s = mean(allData.binnedSpikes);
%         allData.binnedSpikes = bsxfun(@rdivide, bsxfun(@minus, allData.binnedSpikes, mn), s);
%         allData.binnedSpikes(:,s==0) = 0;

        allData.binnedSpikes = bsxfun(@rdivide, allData.binnedSpikes, mn); % normalize by dividing by mean, like assuming poisson
        allData.binnedSpikes(:,mn==0) = 0;
        
    end
    
    if ud.smoothSizeT>0
        fprintf(1, 'smoothing...\n');
        gw = myGaussWin(ud.smoothSizeT, round(1/ud.binSize));
        allData.binnedSpikes = conv2(allData.binnedSpikes, gw, 'same');
    end
end

if exist('traces','var') && ~isempty(traces)
    allData.traces = traces;
else
    allData.traces = [];
end
if ~isempty(auxVid)
    allData.auxVid = auxVid;
else
    allData.auxVid = [];
end

figHandle = figure; 

ud.currentRecTime = startTime;
ud.lastRealTime = now;
ud.rate = 1;
ud.playing = true;
ud.figInitialized = false;
ud.pixel = {[1 1]};
ud.pixelTrace = allData.binnedSpikes(:,ud.pixel{1}(1));
ud.numPanels = 2+numel(allData.auxVid);

ud.cax = [-0.4 0.4];

ud.nColors = 5;
% pixColors = hsv(nColors);
ud.pixColors =  ... % default color order
    [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


if exist('movieSaveFilePath') && ~isempty(movieSaveFilePath)
    WriterObj = VideoWriter(movieSaveFilePath);
    WriterObj.FrameRate=35;
    open(WriterObj);
    ud.WriterObj = WriterObj;
    ud.recording = false;
    set(figHandle, 'Name', 'NOT RECORDING');
else
    ud.WriterObj = [];
    ud.recording = false;
end

set(figHandle, 'UserData', ud);

myTimer = timer(...
    'Period',0.033,...
    'ExecutionMode','fixedRate',...
    'TimerFcn',@(h,eventdata)showNextFrame(h,eventdata, figHandle, allData));

set(figHandle, 'CloseRequestFcn', @(s,c)closeFigure(s,c,myTimer));

set(figHandle, 'KeyPressFcn', @(f,k)movieKeyCallback(f, k));

showNextFrame(1, 1,figHandle, allData)

start(myTimer);

function showNextFrame(h,e,figHandle, allData)

ud = get(figHandle, 'UserData');

nColors = ud.nColors; pixColors = ud.pixColors;

cax = ud.cax;

if ud.playing
    ud.currentRecTime = ud.currentRecTime+(now-ud.lastRealTime)*24*3600*ud.rate;
    if ud.currentRecTime+ud.winSize(2)>max(allData.st)-0.1
        ud.currentRecTime = 0.1-ud.winSize(1); %loop around to the start again;
    end
    ud.lastRealTime = now;
    set(figHandle, 'UserData', ud);
end

if ~ud.figInitialized
%     ax = subtightplot(1,ud.numPanels,1, 0.03, 0.03, 0.03);    
    ax = axes();
    set(ax, 'Position', [0.1 0.05 0.6 0.9]);
    ud.ImageAxisHandle = ax;
    %myIm = imagesc(svdFrameReconstruct(allData.U, allData.V(:, ud.currentFrame)));     
    if ud.binSpikes
        
        inclSamp1 = find(allData.tInds>(ud.winSize(1)+ud.currentRecTime),1);
        inclSamps = inclSamp1:inclSamp1+ud.nBinInWin-1;

        thisT = (0:ud.nBinInWin)*ud.binSize+ud.winSize(1);
        myIm = imagesc(thisT, 1:size(allData.binnedSpikes, 2), allData.binnedSpikes(inclSamps,:)');
    end
        
    set(ax, 'YDir', 'normal', 'YTick', []);
    ud.ImageHandle = myIm;
    
    if ud.normalize
        caxis(cax);
        colormap(colormap_BlueWhiteRed);
    else
        g = gray();
        colormap(g(end:-1:1,:));
        caxis([0 cax(2)]);
    end
    %colormap parula
%     colorbar    
    set(myIm, 'ButtonDownFcn', @(f,k)movieCallbackClick(f, k, allData, figHandle));
    hold on;
    q = plot(ax, ud.pixel{1}(2), ud.pixel{1}(1), 'ko', 'MarkerFaceColor', pixColors(1,:));
    ud.pixMarkerHandles(1) = q;        
    box off 
    drawnow;
    
    % initialize any trace plots here. Use subtightplot and axis off (except
    % the bottom one?)     
    nSP = length(allData.traces)+1;
    currTime = ud.currentRecTime;    
    plotHeight = 0.95/nSP;
    for tInd = 1:nSP-1
        %ax = subtightplot(nSP,ud.numPanels,(tInd-1)*ud.numPanels+2, 0.01, 0.01, 0.01);
        ax = axes(); 
        set(ax, 'HitTest', 'off');
        lowerEdge = 0.01+(tInd-1)*plotHeight; 
        set(ax, 'Position', [0.1 lowerEdge 0.6 plotHeight]);
        ud.traceAxes(tInd) = ax;
        thisT = allData.traces(tInd).t;
        inclT = find(thisT>currTime+ud.winSize(1),1):find(thisT<currTime+ud.winSize(2),1,'last');
        q = plot(thisT(inclT), allData.traces(tInd).v(inclT));
        if isempty(q)
            q = plot(0,0);
        end
        ud.traceHandles{tInd} = q;
        axis off
        xlim(currTime+ud.winSize);
        if isfield(allData.traces(tInd), 'lims') && ~isempty(allData.traces(tInd).lims)
            yl = allData.traces(tInd).lims;
            ylim(yl);
        else
            yl = ylim();
        end
        hold on;
        q = plot([currTime currTime], yl, 'k--');
        ud.traceZeroBars(tInd) = q;
        makepretty;
        
        annotation('textbox', get(ax, 'Position'), 'String', allData.traces(tInd).name, ...
            'EdgeColor', 'none', 'FontSize', 14);
    end
    
    % one more for the selected pixel
    %ax = subtightplot(nSP,ud.numPanels,(nSP-1)*ud.numPanels+2, 0.01, 0.01, 0.01);
    ax = axes();
    lowerEdge = 0.01+(nSP-1)*plotHeight; 
    set(ax, 'Position', [0.1 lowerEdge 0.6 plotHeight]);
    ud.traceAxes(nSP) = ax;
    thisT = allData.tInds;
    inclT = find(thisT>currTime+ud.winSize(1),1):find(thisT<currTime+ud.winSize(2),1,'last');
    q = plot(thisT(inclT), ud.pixelTrace(inclT), 'Color', pixColors(1,:));     
    ud.traceHandles{nSP} = q;
    hold on;
    yl = ylim();
    axis off
    q = plot([currTime currTime], yl, 'k--');    
    ud.traceZeroBars(nSP) = q;
    xlim([currTime-5 currTime+5]);
    makepretty;
    
    
    % any auxVids
    if ~isempty(allData.auxVid)
        for v = 1:numel(allData.auxVid)
            %ax = subtightplot(1,ud.numPanels,v+2, 0.01, 0.01, 0.01);
            ax = axes();
            set(ax, 'Position', [0.75 0.05+(v-1)*0.5 0.2 0.45]);
            allData.auxVid(v).f(ax, currTime, allData.auxVid(v).data);
            title(allData.auxVid(v).name);
            ud.auxAxes(v) = ax;
        end
    end
    
    % anatomical data
    if isfield(allData, 'anatData')
        anatData = allData.anatData;
        axAnat = axes();
        set(axAnat, 'Position', [0.03 0.05 0.05 0.9]);
        
        wfSize = 0.1*ones(size(anatData.coords(:,1)));
        hProbeScatter = scatter(anatData.coords(:,1), anatData.coords(:,2), wfSize);
        hProbeScatter.MarkerFaceColor = 'flat';
        
        hold on;
        for b = 1:size(anatData.borders, 1)
            plot([min(anatData.coords(:,1)) max(anatData.coords(:,1))], ...
                anatData.borders.upperBorder(b)*[1 1], 'k');
            plot([min(anatData.coords(:,1)) max(anatData.coords(:,1))], ...
                anatData.borders.lowerBorder(b)*[1 1], 'k');
            ah = annotation('textbox', 'String', anatData.borders.acronym{b});
            ah.EdgeColor = 'none';
            set(ah, 'Parent', axAnat);
            set(ah, 'Position', [max(anatData.coords(:,1))+10, mean([anatData.borders.lowerBorder(b), anatData.borders.upperBorder(b)]), 0.05, 0.05])
        end
        axis(axAnat, 'off')
        
        allData.hProbeScatter = hProbeScatter;
        allData.axAnat = axAnat;
    end
    
    ud.figInitialized = true;
    set(figHandle, 'UserData', ud);
end

if ud.playing
    if ud.binSpikes
        
        inclSamp1 = find(allData.tInds>(ud.winSize(1)+ud.currentRecTime),1);
        inclSamps = inclSamp1:inclSamp1+ud.nBinInWin-1;
        
        set(ud.ImageHandle, 'CData',  allData.binnedSpikes(inclSamps,:)');
    end
    if length(ud.pixel)>length(ud.pixMarkerHandles)
        ud = get(figHandle, 'UserData');
        for newPix = length(ud.pixMarkerHandles)+1:length(ud.pixel)
            q = plot(ud.ImageAxisHandle, ud.pixel{newPix}(2), ud.pixel{newPix}(1), 'ko', 'MarkerFaceColor', pixColors(mod(newPix-1,nColors)+1,:));
            ud.pixMarkerHandles(newPix) = q;
        end
        set(figHandle, 'UserData', ud);
    end
    
    currTime = ud.currentRecTime;
    set(get(ud.ImageAxisHandle, 'Title'), 'String', sprintf('time %.2f, rate %.2f', currTime, ud.rate));              
    
    nSP = length(allData.traces)+1;
    for n = 1:nSP
        ax = ud.traceAxes(n);        
        set(ud.traceZeroBars(n), 'XData', [currTime currTime]);
        
        if n<nSP
            thisT = allData.traces(n).t;
            inclT = find(thisT>currTime+ud.winSize(1),1):find(thisT<currTime+ud.winSize(2),1,'last');
            set(ud.traceHandles{n}, 'XData', thisT(inclT), 'YData', allData.traces(n).v(inclT));
        elseif n==nSP            
            % last plot is the pixels. It'll be a cell array with multiple
            % pixels
            thisT = allData.tInds;
            inclT = find(thisT>currTime+ud.winSize(1),1):find(thisT<currTime+ud.winSize(2),1,'last');
            
            if length(ud.pixel)>length(ud.traceHandles{n}) 
                ud = get(figHandle, 'UserData');
                % there are new pixels, need to initialize
                for newPix = length(ud.traceHandles{n})+1:length(ud.pixel)
                    q = plot(ax, thisT(inclT), ud.pixelTrace(inclT,newPix), 'Color', pixColors(mod(newPix-1,nColors)+1,:), 'LineWidth', 2.0);
                    ud.traceHandles{n}(newPix) = q;
                end
                set(figHandle, 'UserData', ud);
            end
            
            
            for tr = 1:length(ud.traceHandles{n})
                thisHand = ud.traceHandles{n}(tr);
                set(thisHand, 'XData', thisT(inclT), 'YData', ud.pixelTrace(inclT,tr));
            end
            ylim(ax, cax);
%             mx =  max(ud.pixelTrace(inclT));
%             mn = min(ud.pixelTrace(inclT));
%             if mx>mn
%                 ylim(ud.traceAxes(n), [mn mx]);
%             end
                
        end
        xlim(ax, [currTime+ud.winSize(1) currTime+ud.winSize(2)]);
        
    end
    
    % any auxVids
    if ~isempty(allData.auxVid)
        for v = 1:numel(allData.auxVid)
            ax = ud.auxAxes(v);
            allData.auxVid(v).f(ax, currTime, allData.auxVid(v).data);
        end
    end
    
    drawnow;
    
    
    if ~isempty(ud.WriterObj) && ud.recording
        frame = getframe(figHandle);
        writeVideo(ud.WriterObj,frame);
    end
end






function movieKeyCallback(figHandle, keydata)

ud = get(figHandle, 'UserData');


if ismember(lower(keydata.Key), {'control', 'alt', 'shift'})
    % this happens on the initial press of these keys, so both the Modifier
    % and the Key are one of {'control', 'alt', 'shift'}
    return;
end

ud = get(figHandle, 'UserData');
switch lower(keydata.Key)
    case 'rightarrow'
        ud.currentRecTime = ud.currentRecTime+ud.rate;        
    case 'leftarrow'
        ud.currentRecTime = ud.currentRecTime-ud.rate;
        if ud.currentRecTime<0.1
            ud.currentRecTime = 0.1;
        end
    case 'uparrow'
        ud.rate = ud.rate*1.25;
    case 'downarrow'
        ud.rate = ud.rate/1.25;
    case 'p' % play/pause
        ud.lastRealTime = now; % so that it picks up where it left off when you un-pause
        set(figHandle, 'UserData', ud);
        ud.playing = ~ud.playing;        
    case 'r' %start/stop recording
        ud.recording = ~ud.recording;
        if ~isempty(ud.WriterObj) && ud.recording
            set(figHandle, 'Name', 'RECORDING');
        else
            set(figHandle, 'Name', 'Not recording.');
        end
    case 'c' % clear pixels
        ud.pixel = {ud.pixel{end}};
        ud.pixelTrace = ud.pixelTrace(:,end);
        oldHands = ud.traceHandles{end};
        ud.traceHandles{end} = oldHands(end);
        delete(oldHands(1:end-1));
        oldHands = ud.pixMarkerHandles;
        ud.pixMarkerHandles = oldHands(end);
        delete(oldHands(1:end-1));
        set(ud.pixMarkerHandles, 'MarkerFaceColor', ud.pixColors(1,:));
        set(ud.traceHandles{end}, 'Color', ud.pixColors(1,:));
    case 'hyphen' % scale cax down
        ud.cax = ud.cax/1.25;
        caxis(ud.ImageAxisHandle, ud.cax);
    case 'equal' % scale cax up
        ud.cax = ud.cax*1.25;
        caxis(ud.ImageAxisHandle, ud.cax);
    case 'b'
        ud.currentRecTime = ud.currentRecTime+ud.winSize(1);  
        if ud.currentRecTime<0.1
            ud.currentRecTime = 0.1;
        end
    case 'g' % go to a specified time in the recording
        
        
end
set(figHandle, 'UserData', ud);

function movieCallbackClick(f, keydata, allData, figHandle)

% clickX = keydata.IntersectionPoint(1);
clickY = keydata.IntersectionPoint(2);

clickX = 0;

pixel = round([clickY clickX]);

% fprintf(1, 'new pixel %d, %d\n', pixel(1), pixel(2));

thisPixelTrace = allData.binnedSpikes(:,pixel(1));

ud = get(figHandle, 'UserData');
if keydata.Button == 3
    % new pixel, leave the old one
    ud.pixel{end+1} = pixel;
    ud.pixelTrace(:,end+1) = thisPixelTrace;
elseif keydata.Button == 1
    ud.pixel{end} = pixel;
    ud.pixelTrace(:,end) = thisPixelTrace;
    % update the plotted spot
    set(ud.pixMarkerHandles(end), 'XData', pixel(2), 'YData', pixel(1));
end
set(figHandle, 'UserData', ud);


function closeFigure(s,c,myTimer)
stop(myTimer)
delete(myTimer);

ud = get(s, 'UserData');
if ~isempty(ud.WriterObj)
    close(ud.WriterObj);
end

delete(s);