

function simpleWFplot(wf, coords, nWFtoPlot)

wf = wf(:,20:end-5);
wfAmps = max(wf')-min(wf');
[~,maxCh] = max(wfAmps);
chDists = ((coords(maxCh,1)-coords(:,1)).^2+(coords(maxCh,2)-coords(:,2)).^2).^0.5;
[~,nearestCh] = sort(chDists);
for n = 1:nWFtoPlot
    plot((0:size(wf,2)-1)*0.4+coords(nearestCh(n),1), ...
        wf(nearestCh(n),:)*7+coords(nearestCh(n),2), ...
        'k', 'LineWidth', 2.0); 
    hold on;
end
axis off; 

return;

%% e.g. 

wf = squeeze(sp.waveforms(sp.cids==cid,:,:))';
nWFtoPlot = 16;

simpleWFplot(wf, sp.coords, nWFtoPlot);