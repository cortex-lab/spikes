function readOEwriteDat(oe_path,sync_channel,sync_input)
% function readOEwriteDat(oe_path,sync_channel,sync_input)
%
% Get all .continuous files in folder, break into CH and ADC
% (note: AUX channels are currently ignored, they're for headstage
% accelerometers?)
%
% oe_path - the directory with open-ephys data to convert
% saves CH (ephys.dat) and ADC (sync.dat) files seperately in oe_path
%
% sync_channel - which ADC/TTL channel contains synchronization information 
% NOTE: input TTL channel as 1-indexed like in GUI, not 0-indexed like in
% saved data
%
% sync_input - 'ttl' or 'adc' depending on which input used for sync


%% Set paths and memory
oe_dir = dir([oe_path]);
oe_ch = cellfun(@any, regexp({oe_dir.name},'CH\d.*?.continuous'));

oe_ch_filenames = sort(cellfun(@(x) [oe_path filesep x], ...
    {oe_dir(oe_ch).name},'uni',false));

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
        
        if chunkInd>1 || n>1 % can skip this for the first channel, it was loaded above
            [data, timestamps, info] = load_open_ephys_data_faster(oe_ch_filenames{n});        
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

%% Save synchronization input

switch sync_input
    
    case 'ttl'
        
        % Get events filename
        ttl_file = cellfun(@any, regexp({oe_dir.name},'all_channels.*?.events'));
        ttl_filename = [oe_path filesep oe_dir(ttl_file).name];
        [data, timestamps, info] = load_open_ephys_data_faster(ttl_filename);
        
        % Structure of file: data = channel
        sync_timestamps = timestamps(data == (sync_channel - 1));
        
        % First timestamp is arbitrary with regards to when open ephys was
        % started, so subtract the first timestamp (should be the start of
        % recording timestamp, and be equal to the continuous file 
        % timestamps(1)/sample rate)
        sync_timestamps = sync_timestamps - timestamps(1);
        
        outFile = [oe_path filesep 'sync.mat'];
        save(outFile,'sync_timestamps');
        
    case 'adc'
               
        % Get ADC sync channel filename
        oe_adc = cellfun(@any, strfind({oe_dir.name},['ADC' num2str(sync_channel)]));
        oe_adc_filename = [oe_path filesep oe_dir(oe_adc).name];
        
        chunkSizeSamps = nint16s;
        nChunks = ceil(totalSamps/chunkSizeSamps);
        
        outFile = [oe_path filesep 'sync.dat'];
        fid = fopen(outFile, 'w');
        
        for chunkInd = 1:nChunks
            fprintf(1, 'chunk %d/%d\n', chunkInd, nChunks);
            chunkStart = (chunkInd-1)*chunkSizeSamps+1;
            chunkEnd = min(chunkInd*chunkSizeSamps, totalSamps);
            allData = zeros(nADC, chunkEnd-chunkStart+1, 'int16');
            
            fprintf(1, 'reading ADC %d/%d\n', sync_channel, nADC);
            
            if chunkInd>1 || n>1 % can skip this for the first channel, it was loaded above
                [data, timestamps, info] = load_open_ephys_data_faster(oe_adc_filename);
            end
            
            allData(n,:) = int16(data);
            
            
            fprintf(1, 'writing dat file\n');
            %     outFile = sprintf('%s.dat', filebase, n);
            %     writeDat(outFile, allData)
            fwrite(fid, allData, 'int16');
            fprintf(1, 'done\n')
        end
        
        fprintf(1, 'all chunks finished\n');
        fclose(fid);
        
end
