function readOEwriteDat(oe_path)
% function readOEwriteDat(filebase, nChans)
%
% oe_path is the directory with open-ephys data to convert
% saves CH (ephys.dat) and ADC (sync.dat) files seperately in oe_path

% Get all .continuous files in folder, break into CH and ADC
% (note: AUX channels are currently ignored, they're for headstage
% accelerometers?)#

oe_dir = dir([oe_path filesep '*.continuous']);
oe_ch = cellfun(@any, regexp({oe_dir.name},'CH\d*.continuous'));
oe_adc = cellfun(@any, regexp({oe_dir.name},'ADC\d*.continuous'));

oe_ch_filenames = sort(cellfun(@(x) [oe_path filesep x], ...
    {oe_dir(oe_ch).name},'uni',false));
oe_adc_filenames = sort(cellfun(@(x) [oe_path filesep x], ...
    {oe_dir(oe_adc).name},'uni',false));

% Load in first channel to set memory restrictions
fprintf(1, 'reading chan 1 for initialization\n');
[data, timestamps, info] = load_open_ephys_data_faster(oe_ch_filenames{1});
totalSamps = length(data);

dmem = memory; 
memToLeaveFree = 4 * 2^30; % num of GB to keep free
memToAllocate = dmem.MemAvailableAllArrays - memToLeaveFree;
memToAllocate = max(0, memToAllocate);
nint16s = memToAllocate/2;

%% Save electrophysiology (CH) channels
nCH = length(oe_ch_filenames);
chunkSizeSamps = nint16s/nCH;
nChunks = ceil(totalSamps/chunkSizeSamps);

outFile = [oe_path filesep 'ephys.dat'];
fid = fopen(outFile, 'w');

for chunkInd = 1:nChunks
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunks);
    chunkStart = (chunkInd-1)*chunkSizeSamps+1;
    chunkEnd = min(chunkInd*chunkSizeSamps, totalSamps);
    allData = zeros(nCH, chunkEnd-chunkStart+1, 'int16');

    for n = 1:nCH
        fprintf(1, 'reading CH %d/%d\n', n, nCH);
        thisFile = oe_ch_filenames{n};
        
        if chunkInd>1 || n>1 % can skip this for the first channel, it was loaded above
            [data, timestamps, info] = load_open_ephys_data_faster(thisFile);        
        end

        allData(n,:) = int16(data);

    end

    fprintf(1, 'writing dat file\n');
%     outFile = sprintf('%s.dat', filebase, n);
%     writeDat(outFile, allData)
    fwrite(fid, allData, 'int16');
    fprintf(1, 'done\n')
end

fprintf(1, 'all chunks finished\n');
fclose(fid);

%% Save synchronization input (ADC) channels
nADC = length(oe_adc_filenames);
chunkSizeSamps = nint16s/nADC;
nChunks = ceil(totalSamps/chunkSizeSamps);

outFile = [oe_path filesep 'sync.dat'];
fid = fopen(outFile, 'w');

for chunkInd = 1:nChunks
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunks);
    chunkStart = (chunkInd-1)*chunkSizeSamps+1;
    chunkEnd = min(chunkInd*chunkSizeSamps, totalSamps);
    allData = zeros(nADC, chunkEnd-chunkStart+1, 'int16');

    for n = 1:nADC
        fprintf(1, 'reading ADC %d/%d\n', n, nADC);
        thisFile = oe_adc_filenames{n};
        
        if chunkInd>1 || n>1 % can skip this for the first channel, it was loaded above
            [data, timestamps, info] = load_open_ephys_data_faster(thisFile);        
        end

        allData(n,:) = int16(data);

    end

    fprintf(1, 'writing dat file\n');
%     outFile = sprintf('%s.dat', filebase, n);
%     writeDat(outFile, allData)
    fwrite(fid, allData, 'int16');
    fprintf(1, 'done\n')
end

fprintf(1, 'all chunks finished\n');
fclose(fid);





