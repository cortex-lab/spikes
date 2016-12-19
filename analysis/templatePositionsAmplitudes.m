
function [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)
% function [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)
%
% Compute some basic things about spikes and templates
%
% outputs: 
% - spikeAmps is length nSpikes vector with amplitude in unwhitened space
% of every spike
% - spikeDepths is the position along the probe of every spike (according
% to the position of the template it was extracted with)
% - templateDepths is the position along the probe of every template
% - templateAmps is the amplitude of each template
% - tempsUnW are the unwhitened templates
% - templateDuration is the trough-to-peak time (in samples)
% - waveforms: returns the waveform from the max-amplitude channel
%
% inputs: 
% - temps, the templates (nTemplates x nTimePoints x nChannels)
% - winv, the whitening matrix (nCh x nCh)
% - ycoords, the coordinates of the channels (nCh x 1)
% - spikeTemplates, which template each spike came from (nSpikes x 1)
% - tempScalingAmps, the amount by which the template was scaled to extract
% each spike (nSpikes x 1)

% unwhiten all the templates
tempsUnW = zeros(size(temps));
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
end

% compute the biggest absolute value within each template (obsolete)
% absTemps = abs(tempsUnW);
% tempAmps = max(max(absTemps,[],3),[],2);

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel (but see
% below for true tempAmps)
tempAmpsUnscaled = max(tempChanAmps,[],2);

% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmpsUnscaled*0.3; 
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;

% ... in order to compute the depth as a center of mass
templateDepths = sum(bsxfun(@times,tempChanAmps,ycoords'),2)./sum(tempChanAmps,2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes (templates are zero-indexed)
spikeAmps = tempAmpsUnscaled(spikeTemplates+1).*tempScalingAmps;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
ta = clusterAverage(spikeTemplates+1, spikeAmps);
tids = unique(spikeTemplates);
tempAmps(tids+1) = ta; % because ta only has entries for templates that had at least one spike
tempAmps = tempAmps'; % for consistency, make first dimension template number

% Each spike's depth is the depth of its template
spikeDepths = templateDepths(spikeTemplates+1);

% Get channel with largest amplitude, take that as the waveform
[~,max_site] = max(max(abs(temps),[],2),[],3);
templates_max = nan(size(temps,1),size(temps,2));
for curr_template = 1:size(temps,1)
    templates_max(curr_template,:) = ...
        temps(curr_template,:,max_site(curr_template));
end
waveforms = templates_max;

% Get trough-to-peak time for each template
[~,waveform_trough] = min(templates_max,[],2);
[~,templateDuration] = arrayfun(@(x) ...
    max(templates_max(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(templates_max,1)));





