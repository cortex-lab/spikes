% Inputs: spikeTimes, spikeAmps, spikeYpos - names self explanatory
%         opt - optional, empty by default; 'mark' - will mark detected drifts, 'show' - will generate a different plot, 
%               where only large spikes are used, and the detection of drift locations is demonstrated
function plotDriftmap(spikeTimes, spikeAmps, spikeYpos, opt)
if nargin < 4
  opt = '';
end

if ~strcmpi(opt, 'show')
  nColorBins = 20;
  ampRange = quantile(spikeAmps, [0.1 0.9]);
  colorBins = linspace(ampRange(1), ampRange(2), nColorBins);
  
  colors = gray(nColorBins); colors = colors(end:-1:1, :); % first bin is smalles spikes, starts white
  for b = 1:nColorBins-1
    theseSpikes = spikeAmps>=colorBins(b) & spikeAmps<=colorBins(b+1);
    
    plot(spikeTimes(theseSpikes), spikeYpos(theseSpikes), '.', 'Color', colors(b,:));
    hold on;
  end  
  xlabel('time')
  ylabel('y position')
end
if isempty(opt)
  return
end

for d = 0:800:max(spikeYpos) % break the recording into 800 um segments
  tmp = spikeAmps(spikeYpos >= d & spikeYpos < d+800);
  I = spikeAmps > mean(tmp) + 1.5*std(tmp) & spikeYpos >= d & spikeYpos < d+800; % large spikes in current segment
  driftEvents = detectDriftEvents(spikeTimes(I), spikeYpos(I), strcmpi(opt, 'show'));
  if strcmpi(opt, 'mark') && ~isempty(driftEvents)
    plot(driftEvents(:,1), driftEvents(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
    % text(driftEvents(:,1)+1, driftEvents(:,2), num2str(round(driftEvents(:,3))), 'Color', 'r') % the magnitude of the drift
  end
end
if strcmpi(opt, 'show')
  ylim([min(spikeYpos) max(spikeYpos)]) 
end

% driftEvents will contain a column of times, a column of depths, and a column of drift magnitudes
function driftEvents = detectDriftEvents(spikeTimes, spikeDepths, doPlot)
if nargin < 3
  doPlot = false;
end
driftEvents = [];
if isempty(spikeTimes)
  return % null input ==> nothing to do
end

D = 2; % um
bins = min(spikeDepths)-D:D:max(spikeDepths)+D;
h = histc(spikeDepths, bins);

h = h(1:end-1); % last bin represents the scalar value bins(end), not an interval
bins = bins(1:end - 1) + D/2; % now it's the centre of each interval

if numel(h) < 3
  locs = []; % findpeaks needs an input with >=3 values
else
  [~, locs] = findpeaks(h);
end

if doPlot
  ax(1) = subplot(1, 5, 1); hold on;
  plot(h, bins, 'k')
  box off
  ylabel('y position')
  ax(2) = subplot(1, 5, 2:5); hold on;
  plot(spikeTimes, spikeDepths, '.', 'Color', 0.5*[1 1 1])
  linkaxes(ax, 'y')
  xlim([0 spikeTimes(end)+1])
  xlabel('time')
end

for p = 1:numel(locs)
  if h(locs(p)) < 0.3*spikeTimes(end)
    continue
    % we want the peaks to correspond to some minimal firing rate (otherwise peaks by very few spikes will be considered as well...)
  end
  if doPlot
    subplot(1, 5, 1); hold on;
  end
  
  posBegin = find(h(1:locs(p)) < 0.05*h(locs(p)), 1, 'last');
  if isempty(posBegin)
    posBegin = 1;
  end
  posEnd   = find(h(locs(p):end) < 0.05*h(locs(p)), 1, 'first') + locs(p) - 1;
  if isempty(posEnd)
    posEnd = numel(bins);
  end
  if (p > 1 && posBegin < locs(p-1)) || (p < numel(locs) && posEnd > locs(p+1))
    if doPlot
      plot(h(locs(p)), bins(locs(p)), 'bo')
    end
    continue % no clean enough separation from neighbour peak(s
  elseif doPlot
    plot(h(locs(p)), bins(locs(p)), 'ro')
    plot(xlim, bins(posBegin)*[1 1], '--', 'Color', 0.5*[1 1 1])
    plot(xlim, bins(posEnd)*[1 1], '--', 'Color', 0.5*[1 1 1])
  end
  
  I = spikeDepths > bins(posBegin) & spikeDepths < bins(posEnd);
  
  currentspikeDepths = spikeDepths(I);
  currentspikeTimes  = spikeTimes(I);
  for t = 0:10:spikeTimes(end)-10
    I = currentspikeTimes >= t & currentspikeTimes <= t+10;
    driftSize = bins(locs(p)) - median(currentspikeDepths(I)); 
    if abs(driftSize) > 6 && sum(I) > 10 % 6 um is the hardcoded threshold for drift, and we want at least 10 spikes for the median calculation
      driftEvents(end+1,:) = [t+5, bins(locs(p)), driftSize];
    end
  end
  if doPlot && ~isempty(driftEvents)
    subplot(1, 5, 2:5); hold on;
    plot(driftEvents(:,1), driftEvents(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
    text(driftEvents(:,1)+1, driftEvents(:,2), num2str(round(driftEvents(:,3))), 'Color', 'r') % the magnitude of the drift
  end    
end % loop on peak locations
