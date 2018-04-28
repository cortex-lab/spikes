

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





%% Looking at PSTHs aligned to some event


%% 