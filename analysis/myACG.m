

function [xLin, nLin, xLog, nLog] = myACG(st, axLin, axLog)
% function [xLin, nLin, xLog, nLog] = myACG(st, axLin, axLog)
%
% Computes an autocorrelogram with both linear and log bins
%
% axLin and axLog are optional but if provided will produce a plot in each

binSize = 0.0005;

b = 0.0001:binSize:1; % start at a tenth of a ms to avoid counting itself

[n,xLin] = histdiff(st, st, b);
nLin = n./binSize; 
nLin = nLin(:); xLin = xLin(:);

% now re-bin in a logarithmic way to de-noise at higher lags
q = logspace(-4, 0, 100); % bin edges
xLog = q(1:end-1)+diff(q)/2;

wr = WithinRanges(xLin, [q(1:end-1)' q(2:end)'], 1:numel(q)-1, 'vector');
sums = full(sparse(wr, ones(size(wr)), nLin)); % sum up the acg within each new bin
[v,i] = countUnique(wr); % determine how many of the linear bins went into each log bin
sums(v) = sums(v)./i; % divide by that number so in the end we took a mean of all the linear bins within each log bin
xLog = xLog(ismember(1:numel(xLog), wr)); % some of the log bins didn't contain any linear bins so we drop them
nLog = sums; 
nLog = nLog(ismember(1:numel(nLog), wr));
nLog = nLog(:); xLog = xLog(:);

refractLine = 0.002;    
asymptVal = nLog(end);
    
if ~isempty(axLin)
    plotUntil = 0.1;
    
    b = b(:);
    xx = [b(1); reshape([b(2:end-1) b(2:end-1)]', (numel(b)-2)*2,1); b(end)];
    yy = reshape([nLin nLin]', numel(nLin)*2,1);
    incl = xx<plotUntil;
    plot(axLin, xx(incl), yy(incl), 'k', 'LineWidth', 2.0);
    hold(axLin, 'on');
    plot(axLin, refractLine*[1 1], [0 max(yy(incl))], 'k--');
    plot(axLin, [0 plotUntil], asymptVal*[1 1], 'k--');
    xlim(axLin, [0 plotUntil]);
    if max(yy(incl))>0
        ylim(axLin, [0 max(yy(incl))]);
    end
    set(axLin, 'XTick', [0 0.01 0.1], 'XTickLabel', [0 10 100]);
    set(axLin, 'YTick', []);
    box(axLin, 'off');
end

if ~isempty(axLog)
    bLog = [b(1); xLog(1:end-1)+diff(xLog)/2; q(end)]; %edges
    xx = [bLog(1); reshape([bLog(2:end-1) bLog(2:end-1)]', (numel(bLog)-2)*2,1); bLog(end)];
    yy = reshape([nLog nLog]', numel(nLog)*2,1);
    
    semilogx(axLog, xx, yy, 'k', 'LineWidth', 2.0);
    hold(axLog, 'on');
    semilogx(axLog, refractLine*[1 1], [0 max(yy)], 'k--');
    semilogx(axLog, [min(xx) max(xx)], asymptVal*[1 1], 'k--');
    set(axLog, 'YTick', []);
    set(axLog, 'XTick', [1e-3 1e-2 1e-1 1e0], 'XTickLabel', [1 10 100 1000]);
    if max(yy)>0
        ylim(axLog, [0 max(yy)]);
    end
    xlim([4e-4 max(xx)]);
    box(axLog, 'off');
end