
function [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps)

tempsUnW = zeros(size(temps));
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
end

absTemps = abs(tempsUnW);
% tempChanAmps = squeeze(max(absTemps, [], 2));
tempAmps = max(max(absTemps,[],3),[],2);

tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));
tempAmps = max(tempChanAmps,[],2);

% need to zero-out the potentially-many low values on distant channels
threshVals = tempAmps*0.3; 
tempChanAmps(bsxfun(@lt, tempChanAmps, threshVals)) = 0;

templateYpos = sum(bsxfun(@times,tempChanAmps,ycoords'),2)./sum(tempChanAmps,2);


spikeAmps = tempAmps(spikeTemplates+1).*tempScalingAmps;

spikeDepths = templateYpos(spikeTemplates+1);