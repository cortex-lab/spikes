

function [rmsPerChannel, madPerChannel] = computeRawRMS(ksDir, gainFactor, rawDir)
% function [rmsPerChannel, madPerChannel] = computeRawRMS(ksDir, gain[, rawDir])
%  
% for spikeglx: gainFactor = 0.6/512/gainSetting*1e6;
%
% ksDir is directory of kilosort results
% gainFactor will multiply the raw data
% rawDir is location of raw file, if different from ksDir

pars = loadParamsPy(fullfile(ksDir, 'params.py'));
Fs = pars.sample_rate;
nCh = pars.n_channels_dat;

if nargin>2
    rawFilename = fullfile(rawDir, pars.dat_path);
else
    rawFilename = fullfile(ksDir, pars.dat_path);
end

% read 10 sec of data
% future version, it would be better to memory map it and pull out 10
% one-sec segments (or longer) from various points in the file. But this is
% probably ok, doubt it makes much difference, haven't tested.
fid = fopen(rawFilename, 'r');
rawDat = fread(fid, [nCh Fs*10], 'int16=>double');
fclose(fid);

chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));

% scale by gain
rawDat = rawDat(chanMap+1,:).*gainFactor;

% compute MAD
madPerChannel = mad(rawDat',1)';
rmsPerChannel = rms(rawDat',1)';