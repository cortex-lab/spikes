% Inputs/outputs: mostly self explanatory 
% localizedSpikesOnly (false by default) - if true, only spikes with no discrepancy between depth and site are returned. 
function [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksDir, localizedSpikesOnly)

if nargin < 2
  localizedSpikesOnly = false;
end

clear params
params.loadPCs = true;
sp = loadKSdir(ksDir, params);

if localizedSpikesOnly % go over all templates and check which are not localized (in space)
  localizedTemplates = false(size(sp.temps,1), 1);
  for t = 1:size(sp.temps,1)
    M = max(max(abs(squeeze(sp.temps(t,:,:)))));
    ch = find(max(abs(squeeze(sp.temps(t,:,:)))) > 0.5*M); % the channels where the template has significant weight
    localizedTemplates(t) = max(ch) - min(ch) <= 20; % all channels at most 20 apart?
    
    %assert(sum(sum(squeeze(sp.temps(t,:,:)).^2)) == 1, 'template norm not 1?!')
    %com = sp.ycoords .* max(abs(squeeze(sp.temps(t,:,:)))).^2 / sum(max(abs(squeeze(sp.temps(t,:,:)))).^2)
  end
  localizedTemplates = uint32(find(localizedTemplates) - 1); % the numbers of localized templates (indexing starts from 0)
  i = ismember(sp.spikeTemplates, localizedTemplates);
  sp.st = sp.st(i);
  sp.spikeTemplates = sp.spikeTemplates(i);
  sp.clu = sp.clu(i);
  sp.tempScalingAmps = sp.tempScalingAmps(i);
  sp.pcFeat = sp.pcFeat(i,:,:);
  clear i localizedTemplates M ch t
end

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
spikeFeatYcoords = ycoords(spikeFeatInd+1); % 2D matrix of size #spikes x 12

% center of mass is sum(coords.*features)/sum(features)
spikeDepths = sum(spikeFeatYcoords.*pcFeat.^2,2)./sum(pcFeat.^2,2);


%% for plotting, we need the amplitude of each spike, both so we can draw a
% threshold and exclude the low-amp noisy ones, and so we can color the
% points by the amplitude

[spikeAmps, ~, templateYpos, tempAmps, tempsUnW] = ...
  templatePositionsAmplitudes(temps, winv, ycoords, spikeTemps, tempScalingAmps);

%[~,max_site] = max(max(abs(temps),[],2),[],3); % the maximal site for each template
% one could potentially use the unwhitened templates, but that shouldn't really change the results
[~,max_site] = max(max(abs(tempsUnW),[],2),[],3); % the maximal site for each template
spikeSites = max_site(spikeTemps+1);

if isfield(sp, 'gain') && ~isempty(sp.gain) % could put this field in your params.py
  % spikeAmps = spikeAmps*0.6/512/500*1e6; % convert to uV
  spikeAmps = spikeAmps*sp.gain; % convert to uV
end

if localizedSpikesOnly  % above we already removed non-localized templates, but that on its own is insufficient
  b = regress(spikeSites, spikeDepths); % for IMEC probe adding a constant term kills the regression making the regressors rank deficient
  i = abs(spikeSites - b*spikeDepths) <= 5;
  spikeTimes  = spikeTimes(i);
  spikeAmps   = spikeAmps(i);
  spikeDepths = spikeDepths(i);
  spikeSites  = spikeSites(i);
end
spikeSites = uint16(spikeSites);
