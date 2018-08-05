

function predData = kilosortPredictData(predictSamps, temps, spikeTimes, ...
    spikeTemplates, spikeAmps, invWhiteMat)
% function predData = kilosortPredictData(predictSamps, temps, spikeTimes, ...
%     spikeTemplates, spikeAmps, invWhiteMat)
%
% Given kilosort output, determine what kilosort thought the raw data
% should be based on its model. See cell at the bottom of this file for an
% example usage. 
%
% - predictSamps: [2,1] vector giving start and end sample numbers to be
% predicted
% - temps: [nTemplates, nTimesamples, nChannels] template shapes
% - spikeTimes: [nSpikes, 1] times of all spikes in samples
% - spikeTemplates: [nSpikes, 1] template associated with each spike
% - spikeAmps: [nSpikes, 1] amplitude of each spike (scaling factor of the
% template)
% - invWhiteMat: [nCh, nCh] inverse whitening matrix, to put the templates
% back into real data space

peakShift = 40; % offset of the peak from the start of the waveform
% unfortunately this varies in different versions of kilosort, so you may
% have to try different numbers here. Try (n/2)-1 where n is number of time
% samples in the waveform. 

opsScaleProc = 1; % this is found in rez.ops.scaleproc, and is necessary to
% include here in some versions of kilosort

predictSamps = predictSamps(1):predictSamps(end);
buff = size(temps,2); % width of a spike in samples
nCh = size(temps,3);

% choose the spikes to work with
inclSpikes = spikeTimes>predictSamps(1)-buff/2 & spikeTimes<predictSamps(end)+buff/2;
st = double(spikeTimes(inclSpikes)); % in samples
inclTemps = double(spikeTemplates(inclSpikes))+1;
amplitudes = spikeAmps(inclSpikes);

% un-whiten the waveforms
tempsUnW = zeros(size(temps));
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*invWhiteMat;
end

% add each template shape to the data at the correct time and scaling
predData = zeros(nCh,numel(predictSamps)+buff*4);
for s = 1:sum(inclSpikes)
    
    theseSamps = st(s)+(1:buff)-peakShift-predictSamps(1)+buff*2;
    predData(:,theseSamps) = predData(:,theseSamps) + ...
        squeeze(tempsUnW(inclTemps(s),:,:))' * amplitudes(s);
    
end

% cut buffer and scale
predData = predData(:,buff*2+1:end-buff*2).*opsScaleProc;

return;

%% example usage

ksDir = 'J:\temp';
samps = [0 1000]+60000; 


a = readNPY(fullfile(ksDir, 'amplitudes.npy'));
st = readNPY(fullfile(ksDir, 'spike_times.npy'));
stemps = readNPY(fullfile(ksDir, 'spike_templates.npy'));
temps = readNPY(fullfile(ksDir, 'templates.npy'));
winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));
p = loadParamsPy(fullfile(ksDir, 'params.py'));

predData = kilosortPredictData(samps, temps, st, stemps, a, winv);

figure; 
subplot(1,3,1); 
imagesc(predData);
cax = caxis(); cax = [-1 1]*max(abs(cax))*0.5; caxis(cax);
% colorbar
colormap(colormap_blueblackred); % a divergent colormap
title('predicted data');

fname = fullfile(ksDir, p.dat_path);
d = dir(fname); 
nSamp = d.bytes/2/p.n_channels_dat;
mmf = memmapfile(fname, 'Format', {'int16', [p.n_channels_dat nSamp], 'x'});
realData = mmf.Data.x(:,samps(1):samps(end));
chMap = readNPY(fullfile(ksDir, 'channel_map.npy'))+1;
realData = double(realData(chMap,:));

subplot(1,3,2); 
imagesc(realData)
caxis(cax); 
% colorbar
title('real data');

subplot(1,3,3);
imagesc(realData-predData);
% colorbar
caxis(cax); 
title('real minus prediction');