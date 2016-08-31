

function plotWaveform(wf, xcoords, ycoords, xScale, yScale, thresh, color, varargin)
% wf is nChan x nTimePoints
% color must be 3-element vector

if ~isempty(thresh)
    chanAmps = max(wf,[],2)-min(wf, [], 2);
    maxAmp = max(chanAmps);
    inclChans = chanAmps>maxAmp*thresh;
else
    inclChans = true(size(xcoords));
end

params.LineWidth = 2.0;
params.alpha = 1;
if ~isempty(varargin)
    params = varargin{1};
end
color = [color params.alpha];

t = (0:size(wf,2)-1)/size(wf,2)*xScale;

inclChansNums = find(inclChans);
nIncl = sum(inclChans);
for ch = 1:nIncl
    thisCh = inclChansNums(ch);
    thisWF = wf(thisCh,:) * yScale;
    
    plot(t+xcoords(thisCh), thisWF+ycoords(thisCh), 'k', 'Color', color, 'LineWidth', params.LineWidth);
    hold on;
end
