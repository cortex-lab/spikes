function spikeStruct = loadKSdir(ksDir, varargin)
% spikeStruct = loadKSdir(ksDir) loads kilosorted ephys data. Input
% argument should be a path to a folder containing the output kilosort.
% The output spikeStruct contains the following data:
%
% -dat_path, n_channels_dat, dtype, offset, sample_rate, hp_filtered: these
% are parameter variables associated with the recording. Most variables are
% self-explanatory. 
% -st: vector of all spike event times
% -spikeTemplates: vector of the template identity associated with each
% spike event. [0 indexed]
% -clu: vector of cluster IDs associated with each spike event. This can
% differ from spikeTemplates because spikes from multiple kilosort
% templates can be merged. Or spikes from one template can be split into 
% separate clusters.
% cluster (same value in clu).
% -tempScalingAmps: vector of scaling factors associated with each spike
% event. This value reflects how much the kilosort template for the cluster
% was scaled to match the waveform of each spike event.
% -cids: list of cluster IDs
% -xcoords, ycoords: location of each channel on the probe
% -temps: clusterID x time x channel matrix, containing the kilosort
% template for each cluster.
% -winv: whitening matrix used to ensure all channels have the same
%   variance.
% -pcFeat: PCA analysis on spike statistics
% -pcFeatInd: PCA analysis on spike statistics

if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'excludeNoise')
    params.excludeNoise = true;
end
if ~isfield(params, 'loadPCs')
    params.loadPCs = false;
end

% load spike data

spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
st = double(ss)/spikeStruct.sample_rate;
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed

if exist(fullfile(ksDir, 'spike_clusters.npy'))
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

%Load phy annotation labels (cluster_groups.csv or .tsv)
cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv')) 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
end
if exist(fullfile(ksDir, 'cluster_group.tsv')) 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 
if exist(fullfile(ksDir, 'cluster_KSLabel.tsv')) 
    cgsFile = fullfile(ksDir, 'cluster_KSLabel.tsv');
end
if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);

    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));        
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
        
        
    end
    
else
    clu = spikeTemplates;
    
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end
    

coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);
temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;