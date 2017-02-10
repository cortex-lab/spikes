
function f = plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands)
% function plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands)
%
% Plots LFP power across the probe, depth by frequency, as colormap
% Adds marginals with the power at selected channels, and with the power
% averaged across certain freq bands. 
%
% See companion function lfpBandPower that produces the first two inputs for this.

%%
dispF = F>dispRange(1) & F<=dispRange(2);
nC = size(allPowerEst,1); 

f = figure;
subplot(4,4,[5 6 7 9 10 11 13 14 15]);
imagesc(F(dispF), (0:nC-1)*10, 10*log10(allPowerEst(:,dispF)));
xlim(dispRange);
xlabel('frequency (Hz)');
set(gca, 'YDir', 'normal');
ylabel('depth on probe (µm)');
% h = colorbar;
% h.Label.String = 'power (dB)';
makepretty

ax = subplot(4,4,1:3); hold on;
set(ax, 'ColorOrder', copper(length(marginalChans)));
plot(F(dispF), 10*log10(allPowerEst(marginalChans,dispF)));
ylabel('power (dB)');
set(ax, 'XTick', []);
hleg = legend(array2stringCell(marginalChans*10));
set(hleg, 'Position', [0.7125    0.7607    0.1036    0.2083]);
makepretty;

ax = subplot(4,4,[8 12 16]); hold on;
c = copper(length(freqBands));
c = c(:, [3 2 1]);
set(ax, 'ColorOrder', c);
for q = 1:length(freqBands)
    inclF = F>freqBands{q}(1) & F<=freqBands{q}(2);
    thisPow = mean(10*log10(allPowerEst(:,inclF)),2);
    plot(thisPow, (0:nC-1)*10);
end
set(ax, 'YTick', []);
ylim([0 nC(end)*10])
xlabel('power (dB)');
h = legend(cellfun(@(x)sprintf('%.1f - %.1f Hz', x(1), x(2)), freqBands, 'uni', false));
set(h, 'Position', [ 0.8071    0.7356    0.0980    0.1009]);
makepretty;
