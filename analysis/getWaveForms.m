function wf = getWaveForms(gwfparams)

% % INPUT
% gwfparams.dataDir = '/path/to/data/';                                         % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';                                              % .dat file containing the raw 
% gwfparams.dataType = 'int16'; gwfparams.dataTypeNBytes = 2;                   % Data type and number of bytes per sample in .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                                                           % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];                                                   % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                                                         % Number of waveforms per unit to pull out
% [gwfparams.spikeClusters, gwfparams.spikeTimes] = phy2mat(gwfparams.dataDir); % Spike times of clusters classified as 'good' in Phy
%
% % OUTPUT
% wf.unitIDs        % Cluster IDs (Phy nomenclature)
% wf.spikeTimeKeeps % Which spike times were used for the waveforms
% wf.waveForms      % Individual waveforms
% wf.waveFormsMean  % Average of all waveforms
% 
% % USAGE
% wf = getWaveForms(gwfparams);

% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);             
filenamestruct = dir(fileName);
nSamp = filenamestruct.bytes/(gwfparams.nCh*gwfparams.dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
chMap = readNPY([gwfparams.dataDir 'channel_map.npy'])+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,gwfparams.nCh,wfNSamples);
waveFormsMean = nan(numUnits,gwfparams.nCh,wfNSamples);
for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    curgwfparams.spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(curgwfparams.spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
    end
    waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:)));
    disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
end

% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;

end
