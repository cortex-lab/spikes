

function plotWaveform(wf, xcoords, ycoords, xScale, yScale, thresh, color)
% wf is nChan x nTimePoints

if ~isempty(thresh)
    chanAmps = max(wf,[],2)-min(wf, [], 2);
    maxAmp = max(chanAmps);
    inclChans = chanAmps>maxAmp*thresh;
else
    inclChans = true(size(xcoords));
end

t = (0:size(wf,2)-1)/size(wf,2)*xScale;

inclChansNums = find(inclChans);
nIncl = sum(inclChans);
for ch = 1:nIncl
    thisCh = inclChansNums(ch);
    thisWF = wf(thisCh,:) * yScale;
    
    plot(t+xcoords(thisCh), thisWF+ycoords(thisCh), 'k', 'Color', color, 'LineWidth', 2.0);
    hold on;
end
