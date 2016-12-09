

function [rmsPerChannel, madPerChannel] = computeRawRMS(ksDir, gainFactor)
% function rmsPerChannel = computeRawRMS(ksDir, gain)
% Estimate is based on median absolute deviation, to get back to MAD divide
% by 0.6745. 
% for spikeglx: gainFactor = 0.6/512/gainSetting*1e6;

pars = loadParamsPy(fullfile(ksDir, 'params.py'));
Fs = pars.sample_rate;
nCh = pars.n_channels_dat;
rawFilename = fullfile(ksDir, pars.dat_path);

% read 10 sec of data
fid = fopen(rawFilename, 'r');
rawDat = fread(fid, [nCh Fs*10], 'int16=>double');
fclose(fid);

chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));

% scale by gain
rawDat = rawDat(chanMap+1,:).*gainFactor;

% compute MAD
madPerChannel = mad(rawDat',1)';
rmsPerChannel = rms(rawDat',1)';