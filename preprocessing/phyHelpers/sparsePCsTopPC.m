

function S = sparsePCsTopPC(pcFeat,pcFeatInd, spikeTemplates)
%function S = sparsePCs(pdFeat,pcFeatInd)
% make a sparse matrix version of the pc features
% for the moment, only the top PC

nPCchans = size(pcFeat,3);
nSpikes = size(pcFeat,1);
nTemplates = size(pcFeatInd,1);
nChannels = double(max(pcFeatInd(:))+1);

thisFeat = 1;

rowInds = repmat([1:nSpikes]', [nPCchans,1]);
colInds = zeros(size(rowInds));

for q = 1:nPCchans
    colInds((q-1)*nSpikes+1:q*nSpikes) = pcFeatInd(spikeTemplates+1,q)+1;
end

pcFeatRS = reshape(squeeze(pcFeat(:,thisFeat,:)), [nSpikes*nPCchans 1]);
S = sparse(rowInds, colInds, double(pcFeatRS), nSpikes, nChannels);