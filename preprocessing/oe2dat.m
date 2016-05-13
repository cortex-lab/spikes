function oe2dat(oe_path,save_path,sync_channel,sync_input)
% oe2dat(oe_path,save_path,sync_channel,sync_input)
%
% Converts OE format recorded from open-ephys into flat binary
% Also saves sync channel and parameters/header info
%
% oe_path - path with open-ephys data to convert
% save_path - path to save converted data
%
% sync_channel - which ADC/TTL channel contains synchronization information 
% (NOTE: input TTL channel as 1-indexed like in GUI, not 0-indexed like in
% saved data)
% sync_input - 'ttl' or 'adc' depending on which input used for sync

warning('This isn''t equivalent to kwik2dat yet: doesn''t save all parameters, doesn''t split/filter/CAR data');

%% Make save directory
if ~exist(save_dir,'dir');
    mkdir(savedir)
end

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

outFile = [save_path filesep 'ephys.dat'];
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
        
        outFile = [save_path filesep 'sync.mat'];
        save(outFile,'sync_timestamps');
        
    case 'adc'
               
        % Get ADC sync channel filename
        oe_adc = cellfun(@any, strfind({oe_dir.name},['ADC' num2str(sync_channel)]));
        oe_adc_filename = [oe_path filesep oe_dir(oe_adc).name];
        
        chunkSizeSamps = nint16s;
        nChunks = ceil(totalSamps/chunkSizeSamps);
        
        outFile = [save_path filesep 'sync.dat'];
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

%% Save parameters/header information in separate file

warning('sample rate and gain not included in header yet because I haven''t looked for where they''re stored');
% it's in the 'info' of the function to load OE, I just don't have sample
% data at the moment
params = {'raw_path',['''' oe_path '''']; ...
    'n_channels',num2str(nCH)};
%    'sample_rate',num2str(sample_rate); ...
%    'gain',num2str(ch_gain)};

param_filename = [save_path filesep 'dat_params.txt'];

formatSpec = '%s = %s \r\n';
fid = fopen(param_filename,'w');
for curr_param = 1:size(params,1)
    fprintf(fid,formatSpec,params{curr_param,:});
end
fclose(fid);


