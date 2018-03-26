

function [syncData, detectedFlips, eventTimes] = loadSync(mouseName, thisDate, tag)
%function syncData = loadSync(mouseName, thisDate, tag)
%
% load the sync channel data from a phase3 probe recorded with spikeglx

% rootE = dat.expPath(mouseName, thisDate, 1, 'main', 'master');
% root = fileparts(rootE);
root = getRootDir(mouseName, thisDate);

if nargin>2 && ~isempty(tag)
    dname = ['ephys_' tag];
%     d = fullfile(root, dname, [mouseName '_' thisDate '_' tag '_g0_t0.imec.lf_sync.dat']);
%     lf = fullfile(root, dname, [mouseName '_' thisDate '_' tag '_g0_t0.imec.lf.bin']);
    dDir = dir(fullfile(root, dname, '*.lf_sync.dat'));
    lfDir = dir(fullfile(root, dname, '*.lf.bin'));
else
    dname = 'ephys';
%     d = fullfile(root, dname, [mouseName '_' thisDate '_g0_t0.imec.lf_sync.dat']);
%     lf = fullfile(root, dname, [mouseName '_' thisDate '_g0_t0.imec.lf.bin']);
    dDir = dir(fullfile(root, dname, '*.lf_sync.dat'));
    lfDir = dir(fullfile(root, dname, '*.lf.bin'));
end

if ~isempty(dDir)
    d = fullfile(root, dname, dDir(1).name);
else
    d = [];
end
lf = fullfile(root, dname, lfDir(1).name);

if exist(d)
    fid = fopen(d);
    syncData = fread(fid, [1 Inf], '*int16');
    fclose(fid);
else
    fprintf(1, 'sync data not found at %s\n', d);
    
    if exist(lf)
        fprintf(1, '  extracting sync now...\n')
        syncData = extractSyncChannel(fullfile(root, dname), 385, 385);
    else
        fprintf(1, '  could not find lfp file either at %s\n', lf);
        syncData = [];
    end
end

if nargout>1
    %syncFs = 2500;
    %detectedFlips = schmittTimes((0:length(syncData)-1)/syncFs, syncData, [-1.8 -1.2]);
    
    % get Fs from meta file
    lfmDir = dir(fullfile(root, dname, '*.lf.meta'));
    metaS = readSpikeGLXmeta(fullfile(root, dname, lfmDir(1).name));
    Fs = metaS.sRateHz;
    
    eventTimes = spikeGLXdigitalParse(syncData, Fs);
    detectedFlips = eventTimes{1}{1}; % channel 1, first element (all events - both on and off)
end
