
function [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)
% function [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)
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

% The template amplitude is the amplitude of its largest channel
tempAmps = max(tempChanAmps,[],2);

% need to zero-out the potentially-many low values on distant channels ...
threshVals = tempAmps*0.3; 
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;

% ... in order to compute the depth as a center of mass
templateDepths = sum(bsxfun(@times,tempChanAmps,ycoords'),2)./sum(tempChanAmps,2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes (templates are zero-indexed)
spikeAmps = tempAmps(spikeTemplates+1).*tempScalingAmps;

% Each spike's depth is the depth of its template
spikeDepths = templateDepths(spikeTemplates+1);