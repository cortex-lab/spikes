

function [pdfs, cdfs, ampBins, depthBins] = computeAndPlotWFampHist(ksDir, varargin)
% function [pdfs, cdfs, ampBins, depthBins] = computeAndPlotWFampHist(kilosortDirectory[, ampBins, depthBins])

if ~isempty(varargin)
    ampBins = varargin{1};
    depthBins = varargin{2};
else
    ampBins = [];
    depthBins = [];
end

% load spike data

Fs = 30000; % should read this from params.py ... need to write that function

ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
st = double(ss)/Fs;
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed

if exist(fullfile(ksDir, 'spike_clusters.npy'))
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

if exist(fullfile(ksDir, 'cluster_groups.csv'))
    [cids, cgs] = readClusterGroupsCSV(fullfile(ksDir, 'cluster_groups.csv'));

    noiseClusters = cids(cgs==0);

    st = st(~ismember(clu, noiseClusters));
    spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
    tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
    clu = clu(~ismember(clu, noiseClusters));
    cgs = cgs(~ismember(cids, noiseClusters));
    cids = cids(~ismember(cids, noiseClusters));
end

%
% load forPRBimecP3opt3
% yc = ycoords(connected); xc = xcoords(connected);
coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
yc = coords(:,2); %xc = coords(:,1);
temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

[spikeAmps, spikeDepths, ~, ~, ~] = ...
    templatePositionsAmplitudes(temps, winv, yc, spikeTemplates, tempScalingAmps);

% here, need to scale to uV...

% wf amp plots

if isempty(depthBins)
    depthBins = 0:40:max(yc);
end
if isempty(ampBins)
    ampBins = linspace(0,prctile(spikeAmps,95), 20);
end

recordingDur = max(st);

[pdfs, cdfs] = computeWFampsOverDepth(spikeAmps, spikeDepths, ampBins, depthBins, recordingDur);
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);