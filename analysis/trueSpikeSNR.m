

function snr = trueSpikeSNR(datPars, spikeTimes)

% todo: 
% - dim reduce then reconstruct the waveform with some svd, to avoid
% getting a false boost from overfitting best would be cross-validating
% (i.e. for each spike use only the mean determined from the other
% spikes). Could do this by subtracting the spike off of the mean before
% projecting, actually.
% - use conv to go back through the full raw data trace and get the
% projection, rather than just picking random times. exclude points around
% the detected spikes but otherwise consider the full distribution. 

%% 1. read the actual spikes in to get the mean and projections onto the
% mean

d = dir(datPars.filename);
nSamp = d.bytes/(datPars.nCh*numel(typecast(cast(0, datPars.dataType), 'uint8')));
mmf = memmapfile(datPars.filename, 'Format', {datPars.dataType, [datPars.nCh nSamp], 'x'});

wfNSamps = numel(datPars.wfWin(1):datPars.wfWin(end));
stToUse = spikeTimes(randi(numel(spikeTimes), [datPars.nSpikesToUse 1]));
wfSamps = ceil(stToUse*datPars.Fs)+datPars.wfWin(1);

wfs = zeros(numel(datPars.chanMap), wfNSamps, numel(wfSamps));
for q = 1:numel(wfSamps)
    tmpWf = mmf.Data.x(:,wfSamps(q):wfSamps(q)+wfNSamps-1);
    wfs(:,:,q) = tmpWf(datPars.chanMap,:);
end

mnWF = mean(wfs, 3); 

% for projections, reshape and just multiply
mnWFlin = reshape(mnWF,1,[]);
pOwn = mnWFlin*reshape(wfs, numel(mnWF), []);

if datPars.makePlots
    figure; 
    subplot(1,2,1); 
    imagesc(mnWF); 
    title('mean waveform');
    xlabel('samples'); 
    ylabel('channels');
    subplot(1,2,2);
    imagesc(bsxfun(@rdivide, mnWF, max(abs(mnWF),[],2)+3));
    title('rescaled');
    drawnow;
end

%% 2. read many other random points and get their projections onto the mean

% now pick random times to look at - but NOT random times around existing
% spikes. 
sampsTaken = ceil(spikeTimes*datPars.Fs)+(datPars.wfWin(1):datPars.wfWin(end));
allSamps = 1:nSamp;
sampsAvail = allSamps(~ismember(allSamps, sampsTaken(:))); 
wfSamps = sampsAvail(randi(numel(sampsAvail), [datPars.nSpikesToUse 1])); 
wfs = zeros(numel(datPars.chanMap), wfNSamps, numel(wfSamps));
for q = 1:numel(wfSamps)
    tmpWf = mmf.Data.x(:,wfSamps(q):wfSamps(q)+wfNSamps-1);
    wfs(:,:,q) = tmpWf(datPars.chanMap,:);
end
pOther = mnWFlin*reshape(wfs, numel(mnWF), []);

%% 3. plot histogram/density, compute SNR by gaussian approximation

snr = (mean(pOwn)-mean(pOther))./std(pOther);

if datPars.makePlots
    figure; 
    bins = linspace(min([pOwn pOther]), max([pOwn pOther]), 200);
    [n,x] = hist(pOwn,bins); 
    plot(x,n); hold on; 
    [n,x] = hist(pOther,bins); 
    plot(x,n);
    title(sprintf('snr = %.2f\n', snr))
    xlabel('projection onto mean waveform');
    legend({'spikes in the cluster', 'other'});
end


return;


%% example call
datPars.nCh = 385;
datPars.dataType = 'int16';
datPars.wfWin = [-30 30];
datPars.Fs = 30000;
datPars.makePlots = true;
datPars.nSpikesToUse = 5000;

% serverRoot = '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-11\ephys_PM';
serverRoot = '\\zubjects.cortexlab.net\Subjects\SS096\2018-03-08\ephys_ZO\';
sortingFolder = 'sortingKS2_Nick_Th_8_8';

d = dir(fullfile(serverRoot, '*ap_CAR.bin'));
datPars.filename = fullfile(serverRoot,d.name);

datPars.chanMap = readNPY(fullfile(serverRoot, 'sorting', 'channel_map.npy'))+1;


st = double(readNPY(fullfile(serverRoot, sortingFolder, 'spike_times.npy')))./datPars.Fs;
clu = readNPY(fullfile(serverRoot, sortingFolder, 'spike_clusters.npy'));

spikeTimes = st(clu==43);

snr = trueSpikeSNR(datPars, spikeTimes)
