
function plotWFampCDFs(pdfs, cdfs, ampBins, depthBins)
figure; 

depthX = depthBins(1:end-1)+mean(diff(depthBins))/2;
ampX = ampBins(1:end-1)+mean(diff(ampBins))/2;

subplot(1,2,1); 
imagesc(ampX, depthX, pdfs)
xlabel('spike amplitude (µV)');
ylabel('depth on probe (µm)');
title('pdf');
set(gca, 'YDir', 'normal');
makepretty

subplot(1,2,2); 
imagesc(ampX, depthX, cdfs)
xlabel('spike amplitude (µV)');
ylabel('depth on probe (µm)');
title('inverse cdf');
set(gca, 'YDir', 'normal');
colorbar
colormap(colormap_greyZero_blackred)
caxis([0 20]);
makepretty