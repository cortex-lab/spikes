

function [spikeTimes, spikeAmps, spikeDepths] = ksDriftmap(ksDir)


clear params
params.loadPCs = true;
sp = loadKSdir(ksDir, params);

ycoords = sp.ycoords;
pcFeat = sp.pcFeat;
pcFeat = squeeze(pcFeat(:,1,:)); % take first PC only
pcFeat(pcFeat<0) = 0; % some entries are negative, but we don't really want to push the CoM away from there. 
pcFeatInd = sp.pcFeatInd;
spikeTemps = sp.spikeTemplates;

temps = sp.temps;
winv = sp.winv;
tempScalingAmps = sp.tempScalingAmps;
spikeTimes = sp.st;

%% compute center of mass of these features

% which channels for each spike? 
spikeFeatInd = pcFeatInd(spikeTemps+1,:);

% ycoords of those channels?
spikeFeatYcoords = ycoords(spikeFeatInd+1);

% center of mass is sum(coords.*features)/sum(features)
spikeDepths = sum(spikeFeatYcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);


%% for plotting, we need the amplitude of each spike, both so we can draw a
% threshold and exclude the low-amp noisy ones, and so we can color the
% points by the amplitude

[spikeAmps, ~, templateYpos, tempAmps, tempsUnW] = ...
    templatePositionsAmplitudes(temps, winv, ycoords, spikeTemps, tempScalingAmps);
spikeAmps = spikeAmps*0.6/512/500*1e6; % convert to uV
