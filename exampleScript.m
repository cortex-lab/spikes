

% Example script for some of the functions in the spikes repository. 
%
% These functions make it easy to work with spiking data from kilosort/phy 
% (or, for many functions, from anywhere)
%
% Please make an issue if this function is missing some dependencies.
%
% See https://github.com/cortex-lab/neuropixels/wiki/Other_analysis_methods
% for more explanation. 

%% add the repositories to your path

addpath(genpath('C:\...\github\spikes'))
addpath(genpath('C:\...\github\npy-matlab'))


%% set paths for where to find your data

myKsDir = 'C:\...\data\myKilosortOutputDirectory';

myEventTimes = load('C:\...\data\someEventTimes.mat'); % a vector of times in seconds of some event to align to

%% Loading data from kilosort/phy easily

sp = loadKSdir(myKsDir)

% sp.st are spike times in seconds
% sp.clu are cluster identities
% spikes from clusters labeled "noise" have already been omitted

%% Plotting a driftmap

[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);

%% basic quantification of spiking plot

depthBins = 0:40:3840;
ampBins = 0:30:min(max(spikeAmps),800);
recordingDur = sp.st(end);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);


%% Plotting some basics about LFPs

lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);

lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

chanMap = readNPY(fullfile(myKsDir, 'channel_map.npy'));
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

% plot LFP power
dispRange = [0 100]; % Hz
marginalChans = [10:50:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);

%% Computing some useful details about spikes/neurons (like depths)

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

%% load synchronization data

syncChanIndex = 385;
syncDat = extractSyncChannel(myKsDir, nChansInFile, syncChanIndex);

eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);

% - eventTimes{1} contains the sync events from digital channel 1, as three cells: 
% - eventTimes{1}{1} is the times of all events
% - eventTimes{1}{2} is the times the digital bit went from off to on
% - eventTimes{1}{2} is the times the digital bit went from on to off

% To make a timebase conversion, e.g. between two probes:
% [~,b] = makeCorrection(syncTimesProbe1, syncTimesProbe2, false);

% and to apply it:
% correctedSpikeTimes = applyCorrection(spikeTimesProbe2, b);


%% Looking at PSTHs aligned to some event

% if you now have a vector of relevant event times, called eventTimes (but
% not the cell array as above, just a vector):

window = [-0.3 1]; % look at spike times from 0.3 sec before each event to 1 sec after

% if your events come in different types, like different orientations of a
% visual stimulus, then you can provide those values as "trial groups",
% which will be used to construct a tuning curve. Here we just give a
% vector of all ones. 
trialGroups = ones(size(eventTimes)); 

psthViewer(sp.st, sp.clu, eventTimes, window, trialGroups);

% use left/right arrows to page through the clusters


%% PSTHs across depth

depthBinSize = 80; % in units of the channel coordinates, in this case µm
timeBinSize = 0.01; % seconds
bslWin = [-0.2 -0.05]; % window in which to compute "baseline" rates for normalization
psthType = 'norm'; % show the normalized version
eventName = 'stimulus onset'; % for figure labeling

[timeBins, depthBins, allP, normVals] = psthByDepth(spikeTimes, spikeDepths, ...
    depthBinSize, timeBinSize, eventTimes, window, bslWin);

figure;
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);


%% Loading raw waveforms

% To get the true waveforms of the spikes (not just kilosort's template
% shapes), use the getWaveForms function:

gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==155);

wf = getWaveForms(gwfparams);

figure; 
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;
