

function [rfMap, stats] = sparseNoiseRF(spikeTimes, stimTimes, stimPositions, params)
% function [rfmap, stats] = sparseNoiseRF(spikeTimes, stimTimes, stimPositions, params)
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
    nBins = length([p.countWindow(1):p.binSize:p.countWindow(2)]);
    thisRF = zeros(nX, nY, nBins);
else
    thisRF = zeros(nX,nY);
end

for x = 1:nX
    for y = 1:nY
        theseStims = stimPositions(:,1)==xPos(x) & stimPositions(:,2)==yPos(y);
        [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(spikeTimes, stimTimes(theseStims), p.countWindow, p.binSize);
        if p.useSVD
            thisRF(x,y,:) = psth;
        else
            thisRF(x,y) = mean(spikeCounts);
        end
    end
end

if p.useSVD
    allPSTH = reshape(thisRF, nX*nY, size(thisRF,3));
    bsl = mean(allPSTH(:,1)); % take the first bin as the baseline
    [U,S,V] = svd(allPSTH - bsl,'econ');
    rfMapVec = U(:,1);
    rfMap = reshape(rfMapVec,nX, nY);
%     timeCourse = V(:,1)';
%     Scalar = S(1,1);
%     Model = rfMapVec*timeCourse*Scalar + bsl;
%     Residual = allPSTH - Model;
else
    rfMap = thisRF;
end

if p.fit2Dgauss
end

if p.makePlots
end

stats = [];