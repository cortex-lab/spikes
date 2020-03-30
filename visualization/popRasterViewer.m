
function popRasterViewer(sp, eventData, traces, auxVid, pars)
%
% Viewer for ephys data collected during the task
%
% Displays:
% - Rasters from multiple probes
% -- Can make raster y-axis into cluster identity (sorted by depth) or probe depth
% -- Can make raster colorization be different things, like depth or spike
% amplitude or some feature
%
% - Behavioral traces (lick and wheel)
%
% - Movies (eye and face)
% -- movies loop over the current window
%
% - Event times
%
% Inputs:
% - sp: struct with
%   - st, [nSp 1] spike times
%   - clu, [nSp 1] cluster labels
%   - cids, [nClusters 1] vector of the unique cluster numbers
%   - cgs, [nClusters 1] vector of groups, 1=mua, 2=good, 3=unsorted -
%   optional, only if you want to hide mua
%   - yAxOrderings, [nOrderings 1] struct array with:
%       - name, string
%       - yPos, [nClusters 1] giving a y-axis position for each cluster
%   - colorings, [nColorings 1] struct array with:
%       - name, string
%       - colors, [nClusters 3] colors
% - auxVids: faceVid, eyeVid
% - traces: struct array with
%   - t, timestamps
%   - v, values
%   - name, string
% - eventData: struct array with
%   - color
%   - name
%   - times
%
% Controls:
% - 'y': choose raster y-axis
% - 'c': choose raster colorization
% - left/right arrow: move forward/back in time
% - up/down: make time window bigger/smaller
% - click on raster: draw a vertical line
% - click on event name: toggle it on/off
% - click on trace name: select it - click again: turn off
% - i/j/k/l: move selected trace up/down or scale it
% - '-'/'=': make raster ticks shorter/taller
% - 't': jump to time
% - 'm': hide mua (NOT IMPL)
% - 's'/'f': movies play slower or faster
% - 'p': play/pause movies

% todo: 
% - add anatomical labels
% - option to drop mua
% - add LFP 
% - add screen rendering
% - add ability to do intervals as shaded regions
% - bug with play/pause button
% - try plotting all of events and traces always for speed?


fprintf(1, 'Controls:\n');
fprintf(1, ' - y: choose raster y-axis\n');
fprintf(1, ' - c: choose raster colorization\n');
fprintf(1, ' - left/right arrow: move forward/back in time\n');
fprintf(1, ' - up/down: make time window bigger/smaller\n');
fprintf(1, ' - click on raster: draw a vertical line\n');
fprintf(1, ' - click on event name: toggle it on/off\n');
fprintf(1, ' - click on trace name: select it - click again: turn off\n');
fprintf(1, ' - i/j/k/l: move selected trace up/down or scale it\n');
fprintf(1, ' - -/=: make raster ticks shorter/taller\n');
fprintf(1, ' - t: jump to time\n');
fprintf(1, ' - m: hide mua (NOT IMPL)\n');
fprintf(1, ' - s/f: movies play slower or faster\n');
fprintf(1, ' - p: play/pause movies\n');

for e = 1:length(eventData)
    eventData(e).visible = true;
end
for t = 1:length(traces)
    traces(t).visible = true;
    traces(t).scale = 0.2/(max(traces(t).v)-min(traces(t).v)); % percentage of window
    traces(t).offset = (t-1)/length(traces)+1/(2*length(traces));
end
ud.eventData = eventData;
ud.traces = traces;
ud.auxVid = auxVid;

if isfield(pars, 'startTime'); startTime = pars.startTime; else; startTime = 10; end

params.windowSize = 5; % for timeWin
params.currT = startTime;
params.rasterScale = 1; % height of ticks
params.window = params.windowSize*0.5*[-1 1]+params.currT;
params.binSize = 0.002;
params.playVid = true;
params.posLims = [0 numel(sp.cids)];
params.yAxInd = 1; % which one is currently used
params.colorInd = 1; 
params.selectedTrace = 1;
params.movieT = params.window(1);
params.isPlaying = true;
params.timerPeriod = round(1/30*1000)/1000;
params.movieRate = 1; % 1 = realtime
params.lastRealTime = now;
params.muaHidden = false;

% re-scale all the y-orderings to be the same so things don't re-scale when
% the ordering covers a different range
for y = 1:numel(sp.yAxOrderings)
    yp = sp.yAxOrderings(y).yPos;
    ylv = sp.yAxOrderings(y).yLabelVal;
    mn = min(yp); mx = max(yp);
%     sp.yAxOrderings(y).yLabels = [mn mx];
    sp.yAxOrderings(y).yPos = (yp-mn)/mx*diff(params.posLims)+params.posLims(1);
    sp.yAxOrderings(y).yLabelVal = (ylv-mn)/mx*diff(params.posLims)+params.posLims(1);
    
    
end
ud.sp = sp;
ud.params = params;

f = figure;

ud.f = f;
hands = createPlots(ud);
ud.hands = hands;
updatePlots(ud);
recolor(ud);

set(f, 'UserData', ud);
set(f, 'KeyPressFcn', @(f,k)keyboardPressCallback(f, k));

myTimer = timer(...
    'Period',params.timerPeriod,...
    'ExecutionMode','fixedRate',...
    'TimerFcn',@(h,eventdata)updateMovies(f));

set(f, 'CloseRequestFcn', @(s,~)closeFigure(s, myTimer));

start(myTimer);

function h = createPlots(ud)
sp = ud.sp; st = sp.st; clu = sp.clu;
ed = ud.eventData;
tr = ud.traces;
p = ud.params;

h = [];
h.rasterAx = axes('Position', [0.03 0.03 0.7 0.96]);
set(h.rasterAx, 'Color', 'k');
set(h.rasterAx, 'ButtonDownFcn', @(~,k)rastClick(k, ud.f));
hold(h.rasterAx, 'on');

for c = 1:length(sp.cids)
    thisH = plot(0,0,'w-', 'LineWidth', 1.5);        
    
%     varX = ['rasterData' num2str(c) 'X']; varY = ['rasterData' num2str(c) 'Y'];
%     h.(varX) = 0; h.(varY) = 0; 
%     set(thisH, 'XDataSource', ['h.' varX], 'YDataSource', ['h.' varY]); 
    
    set(thisH, 'HitTest', 'off');
    h.rasterHands(c) = thisH;
end
box(h.rasterAx, 'off');

for e = 1:length(ed)    
    thisH = plot(0,0,'-',ed(e).spec{:});
    set(thisH, 'HitTest', 'off');
    h.eventHands(e) = thisH;
end

for t = 1:length(tr)
    thisH = plot(0,0,'-', 'Color', tr(t).color);
    set(thisH, 'HitTest', 'off');
    h.traceHands(t) = thisH;
end

h.userLine = plot([0 0], p.posLims, 'w');

h.ctrlPanel = uix.Panel('Parent', ud.f, 'Title', 'Controls', 'Position', [0.75 0.55 0.24 0.44]);
vb = uix.VBox('Parent', h.ctrlPanel);
hb = uix.HBox('Parent', vb);
h.ColorDisplay = uicontrol('Style', 'text', 'HorizontalAlignment', 'left', ...
    'String', '','Parent', hb);
h.yPosDisplay = uicontrol('Style', 'text', 'HorizontalAlignment', 'left', ...
    'String', '','Parent', hb);
hb2 = uix.HBox('Parent', vb);
h.evtsPanel = uix.Panel('Parent', hb2, 'Title', 'Events');
h.trPanel = uix.Panel('Parent', hb2, 'Title', 'Traces');
h.intPanel = uix.Panel('Parent', hb2, 'Title', 'Intervals');
vbe = uix.VBox('Parent', h.evtsPanel);
vbt = uix.VBox('Parent', h.trPanel);
vbi = uix.VBox('Parent', h.intPanel);

for e = 1:length(ed)
    %thisH = annotation('textbox', [0.75 0.8-length(tr)*0.05-e*0.05 0.1 0.1], 'String', ['[E] ' ed(e).name]);
    thisH = uicontrol('Style', 'pushbutton', 'HorizontalAlignment', 'left', ...
    'String', ed(e).name,'Parent', vbe);
    thisH.BackgroundColor = 'k';
    thisH.FontWeight = 'bold';
    ii = find(strcmp(ed(e).spec, 'Color'));
    thisH.ForegroundColor = ed(e).spec{ii+1};
    set(thisH, 'Callback', @(~,~)toggleEvent(ud.f, e))
    set(thisH, 'KeyPressFcn', @(a,k)keyboardPressCallback(ud.f,k));
    h.evLeg(e) = thisH;
end

for t = 1:length(tr)
    %thisH = annotation('textbox', [0.75 0.8-t*0.05 0.1 0.1], 'String', ['[T] ' tr(t).name]);
    thisH = uicontrol('Style', 'pushbutton', 'HorizontalAlignment', 'left', ...
    'String', tr(t).name,'Parent', vbt);
    thisH.BackgroundColor = 'k';
    thisH.ForegroundColor = tr(t).color;
    set(thisH, 'Callback', @(~,~)selectTrace(ud.f, t))
    set(thisH, 'KeyPressFcn', @(a,k)keyboardPressCallback(ud.f,k));
    h.trLeg(t) = thisH;
end
vb.Heights = [25 -1];



% h.ColorDisplay = annotation('textbox', [0.75 0.9 0.1 0.1], 'String', '', 'EdgeColor', 'none');
% h.yPosDisplay = annotation('textbox', [0.75 0.85 0.1 0.1], 'String', '', 'EdgeColor', 'none');
% 
% for e = 1:length(ed)
%     thisH = annotation('textbox', [0.75 0.8-length(tr)*0.05-e*0.05 0.1 0.1], 'String', ['[E] ' ed(e).name]);
%     thisH.BackgroundColor = 'k';
%     thisH.FontWeight = 'bold';
%     ii = find(strcmp(ed(e).spec, 'Color'));
%     thisH.Color = ed(e).spec{ii+1};
%     set(thisH, 'ButtonDownFcn', @(~,~)toggleEvent(ud.f, e))
%     h.evLeg(e) = thisH;
% end
% 
% for t = 1:length(tr)
%     thisH = annotation('textbox', [0.75 0.8-t*0.05 0.1 0.1], 'String', ['[T] ' tr(t).name]);
%     thisH.BackgroundColor = 'k';
%     thisH.Color = tr(t).color;
%     set(thisH, 'ButtonDownFcn', @(~,~)selectTrace(ud.f, t))
%     h.trLeg(t) = thisH;
% end



if ~isempty(ud.auxVid)
    for v = 1:numel(ud.auxVid)
        ax = axes('Position', [0.75 0.05+(v-1)*0.25 0.23 0.25]);
        ud.auxVid(v).f(ax, p.movieT, ud.auxVid(v).data);
        h.auxAxes(v) = ax;
    end
    h.movieTimeLine = plot(h.rasterAx,[0 0], p.posLims, 'w:');
else
    h.movieTimeLine = [];
end


function updateMovies(f)
ud = get(f, 'UserData');
p = ud.params;
if p.isPlaying
    newT = p.movieT+(now-p.lastRealTime)*24*3600*p.movieRate;
    if newT>p.window(2)
        newT = p.window(1);
    end
    p.lastRealTime = now;
    p.movieT = newT;
    ud.params = p;
    set(f, 'UserData', ud);
    if ~isempty(ud.auxVid)
        for v = 1:numel(ud.auxVid)
            ax = ud.hands.auxAxes(v);
            ud.auxVid(v).f(ax, p.movieT, ud.auxVid(v).data);
        end
    end
    if ~isempty(ud.hands.movieTimeLine)
        set(ud.hands.movieTimeLine, 'XData', p.movieT*[1 1]);
        title(ud.hands.auxAxes(v), sprintf('rate = %.2f', p.movieRate));
    end
    drawnow;
end

function updatePlots(ud)

sp = ud.sp;
p = ud.params;
h = ud.hands;
ed = ud.eventData;
tr = ud.traces;
set(h.rasterAx, 'XLim', p.window, 'YLim', p.posLims)

incl = sp.st>p.window(1) & sp.st<p.window(2);
st = sp.st(incl); clu = sp.clu(incl);
for c = 1:length(sp.cids)
    [rasterX,yy] = rasterize(st(clu==sp.cids(c)));
    rasterY = yy*p.rasterScale+sp.yAxOrderings(p.yAxInd).yPos(c); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
    set(h.rasterHands(c), 'XData', rasterX, 'YData', rasterY);
    
%     varX = ['rasterData' num2str(c) 'X']; varY = ['rasterData' num2str(c) 'Y'];
%     h.(varX) = rasterX; h.(varX) = rasterY; 
end
% refreshdata(h.rasterAx, 'caller');

for e = 1:length(ed)
%     if ed(e).visible
        incl = ed(e).times>p.window(1) & ed(e).times<p.window(2);
        if sum(incl)>0
            [xx,yy] = rasterize(ed(e).times(incl));    
            set(h.eventHands(e), 'XData', xx, 'YData', yy*p.posLims(2));
        else
            set(h.eventHands(e), 'XData', [], 'YData', []);
        end
%     else
%         set(h.eventHands(e), 'XData', [], 'YData', []);
%     end
end

for t = 1:length(tr)
%     if tr(t).visible         
        incl = tr(t).t>p.window(1) & tr(t).t<p.window(2);
        set(h.traceHands(t), 'XData', tr(t).t(incl), ...
            'YData', tr(t).v(incl)*tr(t).scale*diff(p.posLims)+tr(t).offset*diff(p.posLims));
%     else
%         set(h.traceHands(t), 'XData', [], 'YData', []);
%     end
end

set(h.yPosDisplay, 'String', sprintf('ordered by %s', sp.yAxOrderings(p.yAxInd).name));
set(h.rasterAx, 'YTickLabel', sp.yAxOrderings(p.yAxInd).yLabel, 'YTick', sp.yAxOrderings(p.yAxInd).yLabelVal); 


function recolor(ud)

sp = ud.sp;
p = ud.params;
h = ud.hands;

cmap = sp.colorings(p.colorInd).colors;
if p.muaHidden && isfield(sp, 'cgs')
    set(h.ColorDisplay, 'String', sprintf('colored by %s (mua off)', sp.colorings(p.colorInd).name));
    cmap(:,4) = 1;
    cmap(sp.cgs<2,4) = 0; % full transparency for anything not Good or Unsorted
else
    set(h.ColorDisplay, 'String', sprintf('colored by %s', sp.colorings(p.colorInd).name));
end
    
for c = 1:length(sp.cids)     
    set(h.rasterHands(c), 'Color', cmap(c,:));
end

function rastClick(keydata, f)
clickX = keydata.IntersectionPoint(1);
ud = get(f, 'UserData');
set(ud.hands.userLine, 'XData', clickX*[1 1]);

function selectTrace(f, t)
ud = get(f, 'UserData');
p = ud.params;
tr = ud.traces;
h = ud.hands;
if p.selectedTrace == t    
%     tr(t).visible = false;
    set(h.traceHands(t), 'Visible', 'off');
    ud.hands.trLeg(t).FontWeight = 'normal';   
    p.selectedTrace = [];
else
    if ~isempty(p.selectedTrace)
        ud.hands.trLeg(p.selectedTrace).FontWeight = 'normal';    
    end
    ud.hands.trLeg(t).FontWeight = 'bold';
    p.selectedTrace = t; 
%     tr(t).visible = true;
    set(h.traceHands(t), 'Visible', 'on');
end
ud.traces = tr;
ud.params = p;
set(f, 'UserData', ud);
% updatePlots(ud);

function toggleEvent(f, e)
ud = get(f, 'UserData');
% p = ud.params;
% ed = ud.eventData;
h = ud.hands;
% if ed(e).visible  
if strcmp(get(h.eventHands(e), 'Visible'), 'on')  
%     ed(e).visible = false;
%     set(h.eventHands(e), 'XData', [], 'YData', []);
    set(h.eventHands(e), 'Visible', 'off');
    ud.hands.evLeg(e).FontWeight = 'normal';   
else
%     ed(e).visible = true;
    set(h.eventHands(e), 'Visible', 'on');
    ud.hands.evLeg(e).FontWeight = 'bold';  
%     ud.eventData = ed;
%     updatePlots(ud);
end
% ud.eventData = ed;
% ud.params = p;
% set(f, 'UserData', ud);

function keyboardPressCallback(f, keydata)

ud = get(f, 'UserData');
p = ud.params;
tr = ud.traces;
switch keydata.Key
    case 'rightarrow'
        p.currT = p.currT+p.windowSize/4;
        if p.currT>ud.sp.st(end); p.currT = 0; end
        p = updateWindow(p);        
        if p.movieT<p.window(1)
            p.movieT = p.window(1);
        end
    case 'leftarrow' 
        p.currT = p.currT-p.windowSize/4;
        if p.currT<0; p.currT = ud.sp.st(end); end
        p = updateWindow(p);    
        if p.movieT>p.window(2)
            p.movieT = p.window(1);
        end
    case 'uparrow'
        p.windowSize = p.windowSize*0.75;
        p = updateWindow(p);
        if p.movieT<p.window(1) || p.movieT>p.window(2)
            p.movieT = p.window(1);
        end
    case 'downarrow'
        p.windowSize = p.windowSize/0.75;
        p = updateWindow(p);
        if p.movieT<p.window(1) || p.movieT>p.window(2)
            p.movieT = p.window(1);
        end
    case 'y'
        p.yAxInd = p.yAxInd+1;
        if p.yAxInd>numel(ud.sp.yAxOrderings); p.yAxInd = 1; end
    case 'c'
        p.colorInd = p.colorInd+1;
        if p.colorInd>numel(ud.sp.colorings); p.colorInd = 1; end
        isP = p.isPlaying;
        p.isPlaying = false; ud.params = p; set(f, 'UserData', ud)
        recolor(ud);
        p.isPlaying = isP; ud.params = p; set(f, 'UserData', ud);

    case 'hyphen'
        p.rasterScale = p.rasterScale/1.5;
    case 'equal'
        p.rasterScale = p.rasterScale*1.5;
    case 't'
        newT = inputdlg('jump to time?');
        newT = str2double(newT{1});
        if ~isempty(newT)
            if newT<0
                p.currT = 0; 
            elseif newT>ud.sp.st(end)
                p.currT = ud.sp.st(end);
            else
                p.currT = newT;
            end
        end
        p = updateWindow(p);
    case 'i'
        tr(p.selectedTrace).offset = tr(p.selectedTrace).offset + 1/25;
    case 'k'
        tr(p.selectedTrace).offset = tr(p.selectedTrace).offset - 1/25;
    case 'l'
        tr(p.selectedTrace).scale = tr(p.selectedTrace).scale*5/4;
    case 'j'
        tr(p.selectedTrace).scale = tr(p.selectedTrace).scale*4/5;
    case 'p' 
        if p.isPlaying
            p.isPlaying = false;
        else
            p.isPlaying = true;
            p.lastRealTime = now;
        end
    case 's'
        p.movieRate = p.movieRate*4/5;
    case 'f'
        p.movieRate = p.movieRate*5/4;
    case 'm'
        p.muaHidden = ~p.muaHidden;
        ud.params = p;
        recolor(ud);
end
ud.traces = tr;

isP = p.isPlaying; 
p.isPlaying = false; 
ud.params = p;
set(f, 'UserData', ud)
updatePlots(ud);

p.isPlaying = isP; 
ud.params = p;
set(f, 'UserData', ud);

function p = updateWindow(p)
p.window = p.windowSize*0.5*[-1 1]+p.currT;

function closeFigure(f, myTimer)
stop(myTimer);
delete(myTimer);
delete(f);
