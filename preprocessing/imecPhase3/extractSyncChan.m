
function syncData = extractSyncChan(lfpFilename)
% syncData = function extractSyncChan(lfpFilename)
% 
% extracts the synchronization channel from a neuropixels phase3 data file
% assumes it is channel 385 out of 385, and that data type is int16

[path, fname, ext] = fileparts(lfpFilename);
syncOut = fullfile(path, [fname '_sync' ext]);
syncNchans = 385;

d = dir(lfpFilename);
nSamps = d.bytes/2/syncNchans;
mmf = memmapfile(lfpFilename, 'Format', {'int16', [syncNchans nSamps], 'x'});

syncData = mmf.Data.x(end,:);

fid = fopen(syncOut, 'w');
fwrite(fid, syncData, 'int16');
fclose(fid);