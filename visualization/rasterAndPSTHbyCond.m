

function rasterAndPSTHbyCond(st, evTimes, evGroups, win, sm, colors, axR, axP)
% sm is smoothing window size in sec
% axR and axP are axes for rasters and psth's

%% - compute psths
binSize = 0.0005;
[~,ii] = sort(evTimes);
[~,revii] = sort(ii);
[~, bins, ~, ~, ~, ba] = psthAndBA(...
    st, evTimes, win+[-0.1 0.1], binSize);
ba = ba(revii,:)./binSize;


%% - plot rasters

grp = unique(evGroups);
nEv = numel(grp);

for e = 1:nEv
    [ii,jj] = find(ba(evGroups==grp(e),:));
    ii = ii+sum(ismember(evGroups, grp(1:e-1)));
    
    [rx, ry] = rasterize(bins(jj));
    ry(1:3:end) = ry(1:3:end)+ii'-1;ry(2:3:end) = ry(2:3:end)+ii'-1;
    
    plot(axR, rx, ry, 'Color', colors(e,:), 'LineWidth', 1.5); 
    hold(axR, 'on');
    
    
    plot(axR, min(win)*[1 1]+0.001, sum(ismember(evGroups, grp(1:e-1)))+[0 sum(evGroups==grp(e))], 'Color', colors(e,:), 'LineWidth', 3); 
end

ylabel(axR, 'trials')

if 0>win(1) && 0<win(end)
    plot(axR, [0 0], [0 numel(evGroups)], 'k', 'LineWidth', 2.0);
end
xlim(axR, win); 
ylim(axR, [0 numel(evGroups)]); 

box(axR, 'off'); set(axR, 'XTick', []);
% axis off;


%% psth for two kinds of responses

mgw = myGaussWin(sm, 1/binSize); mgw(1:round(numel(mgw)/2)-1) = 0; mgw = mgw./sum(mgw);

clear h

baSm = conv2(1,mgw, ba, 'same');
for e = 1:nEv
        
    mn = mean(baSm(evGroups==grp(e),:)); 
    stderr = std(baSm(evGroups==grp(e),:))./sqrt(sum(evGroups==grp(e)));
    h(e) = plotWithErr(bins, mn, stderr, colors(e,:));
%     
%     p = conv(mean(ba(evGroups==grp(e),:)), mgw, 'same')./binSize;
%     
%     h(e) = plot(axP, bins, p, 'Color', colors(e,:), 'LineWidth', 2.0); 
    hold(axP, 'on');
    
end

plot(axP, [0 0], ylim(), 'r');

% legend(h, {'reward', 'neg feedback', 'go cue'}, 'box', 'off',...
%     'Location', 'NorthOutside');

xlabel('time from event (s)');
ylabel('firing rate (sp/s)');
xlim(win); box off;