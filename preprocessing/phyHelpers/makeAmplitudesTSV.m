

function makeAmplitudesTSV(ksDir)
% function makeAmplitudesTSV(ksDir)
% makes a tab-separated value file with the amplitudes of each template in
% the given kilsort results directory

s = loadKSdir(ksDir);

tids = unique(s.spikeTemplates);

[~, ~, ~, tempAmps, ~, ~, ~] = templatePositionsAmplitudes(...
    s.temps, s.winv, s.ycoords, s.spikeTemplates, s.tempScalingAmps);

if isfield(s, 'gain') % convert to µV if available
    tempAmps = tempAmps*s.gain;
end

writePhyTSV(ksDir, 'amp', tids, tempAmps(tids+1));

