

function plotDriftmap(spikeTimes, spikeAmps, spikeYpos)

nColorBins = 20;
ampRange = quantile(spikeAmps, [0.1 0.9]);
colorBins = linspace(ampRange(1), ampRange(2), nColorBins);

colors = gray(nColorBins); colors = colors(end:-1:1, :); % first bin is smalles spikes, starts white
for b = 1:nColorBins-1
    theseSpikes = spikeAmps>=colorBins(b) & spikeAmps<=colorBins(b+1); 
    
    plot(spikeTimes(theseSpikes), spikeYpos(theseSpikes), '.', 'Color', colors(b,:));
    hold on;
    
end


xlabel('time');
ylabel('y position');