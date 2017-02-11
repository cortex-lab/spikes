
function plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType, ax)
% function plotPSTHbyDepth(timeBins, depthBins, allP, eventName, ax)
%
% see matching function "psthByDepth" that generates input for this
%
% psthType is 'norm' if you have normalized units, otherwise 

nD = length(depthBins)-1;

if nargin<6 || isempty(ax)    
    ax = gca;
end

imagesc(timeBins, depthBins(1:end-1), allP);
set(gca, 'YDir', 'normal');
hold on;

if strcmp(psthType, 'norm')
    plot(zeros(1, nD), depthBins(1:end-1), 'k--', 'LineWidth', 2.0)
    colormap(ax, colormap_BlueWhiteRed)
    caxis([-10 10]);
else
    plot(zeros(1, nD), depthBins(1:end-1), 'w--', 'LineWidth', 2.0)
end
    
xlabel(['time from ' eventName ' (sec)']);
xlim([min(timeBins) max(timeBins)]);
ylabel('depth on electrode array (µm)')
box off
h = colorbar;
if strcmp(psthType, 'norm')
    h.Label.String = 'Firing rate z-score';
else
    h.Label.String = 'Firing rate (spikes/sec)';
end
makepretty;

