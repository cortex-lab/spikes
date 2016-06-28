

function [rfMap, stats] = sparseNoiseRF(spikeTimes, stimTimes, stimPositions, params)
% function [rfmap, stats] = sparseNoiseRF(spikeTimes, stimTimes, stimPositions, params)
%
% Assumes that stimPositions describe a rectangle and are evenly spaced.
%
% - stimPositions is Nx2
%
% params can have: 
% - makePlots - logical
% - useSVD - logical, whether to compute the RF by using SVD on the PSTHs for all
% stimulus responses. If not, will just count spikes in a window. 
% - countWindow - 1x2, start and end time relative to stimulus onset to
% consider
% - binSize - 1x1, to use for making rasters
% - fit2Dgauss - logical, try to fit a 2D gaussian

% default parameters
p.makePlots = true; 
p.useSVD = true;
p.countWindow = [0 0.2];
p.binSize = 0.01;
p.fit2Dgauss = false;

fn = fieldnames(p);
for f = 1:length(fn)
    if isfield(params, fn{f}) && ~isempty(params.(fn{f}))
        p.(fn{f}) = params.(fn{f});
    end
end

xPos = unique(stimPositions(:,1)); nX = length(xPos);
yPos = unique(stimPositions(:,2)); nY = length(yPos);

if p.useSVD
    timeBins = [p.countWindow(1):p.binSize:p.countWindow(2)];
    timeBins = timeBins(1:end-1)+p.binSize/2;
    nBins = length(timeBins);
    thisRF = zeros(nX, nY, nBins);
else
    thisRF = zeros(nX,nY);
end

% for x = 1:nX
%     for y = 1:nY
%         theseStims = stimPositions(:,1)==xPos(x) & stimPositions(:,2)==yPos(y);
% %         [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(spikeTimes, stimTimes(theseStims), p.countWindow, p.binSize);
%         [psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(spikeTimes, stimTimes(theseStims), p.countWindow, p.binSize);
%             
%         if p.useSVD
%             thisRF(x,y,:) = psth;
%         else
%             thisRF(x,y) = mean(spikeCounts);
%         end
%     end
% end
[psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(spikeTimes, stimTimes, p.countWindow, p.binSize);

for x = 1:nX
    for y = 1:nY
        theseStims = stimPositions(:,1)==xPos(x) & stimPositions(:,2)==yPos(y);
        if p.useSVD
            thisRF(x,y,:) = mean(ba(theseStims,:));
        else
            thisRF(x,y) = mean(spikeCounts(theseStims));
        end
    end
end

if p.useSVD
    allPSTH = reshape(thisRF, nX*nY, size(thisRF,3));
    bsl = mean(allPSTH(:,1)); % take the first bin as the baseline
    [U,S,V] = svd(allPSTH - bsl,'econ');
    rfMapVec = U(:,1);
    rfMap = reshape(rfMapVec,nX, nY);
    timeCourse = V(:,1)';
    Scalar = S(1,1);
    Model = rfMapVec*timeCourse*Scalar + bsl;
    Residual = allPSTH - Model;
else
    rfMap = thisRF;
end

if p.fit2Dgauss
    
    upSampMap = interpn(rfMap,[5 5], 'cubic'); 
    xcoords = (0:size(upSampMap,1)-1)/(max(xPos(:))-min(xPos(:)))+min(xPos(:));
    ycoords = (0:size(upSampMap,2)-1)/(max(yPos(:))-min(yPos(:)))+min(yPos(:));
    figure;
    [x, fitResp] = fit2dGaussRF(ycoords, xcoords, upSampMap, 1);
    x
end

if p.makePlots
    
    
    figure; 
    
    
    if p.useSVD
        nCol = 2;
        
        subplot(3,nCol,2);
        plot(timeBins, timeCourse);
        xlim([timeBins(1) timeBins(end)]);
        title('response time course')
        
        subplot(3,nCol,4);
        imagesc(allPSTH');
        title('all PSTHs, ordered');
        
        subplot(3,nCol,6); 
        imagesc(Residual');
    else
        nCol = 1;
    end
    
    subplot(3,nCol,1);
    imagesc(rfMap);
    title('map');
    
    if p.fit2Dgauss
        subplot(3,nCol,(nCol-1)*2+1)
        imagesc(upSampMap);
        title('up-sample map');
        
    end
    
    
    
end

stats = [];