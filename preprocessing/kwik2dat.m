function kwik2dat(kwik_path,save_path,sync_channels,sync_input,local_copy,subtract_light)
% kwik2dat(kwik_path,save_path,sync_channel,sync_input)
% IN PROGRESS NOTES: started to set up to load in chunks (by channel for
% filtering/downsampling LFP, but then this couldn't do CAR...) didn't
% finish
%
% Converts kwik format recorded from open-ephys into flat binary
% Splits signal in two: one 500 Hz low-pass LFP(lfp.dat), 
% one cross-channel median subtracted spike (spikes.dat)
% Also saves sync channel and parameters/header info
%
% kwk_path - path with kwik files
% save_path - path to save the converted data
%
% sync_channels - TTL channels with synchronization input to save
% (NOTE: input TTL channel as 1-indexed like in GUI, not 0-indexed like in
% saved data)
% sync_input - 'ttl' or 'adc', depending on digital/analog in of sync
% local_copy - true/false, copies data to temporary local folder for faster reading


% Kwik files from open ephys: 
% .raw.kwd: continuous data from all channels
% .kwe: event data

%% Make save directory
if ~exist(save_path,'dir');
    mkdir(save_path)
end

%% Get filenames for each experiment in path
kwd_dir = dir([kwik_path filesep '*.kwd']);
kwd_filename = [kwik_path filesep kwd_dir(1).name];

kwe_dir = dir([kwik_path filesep '*.kwe']);
kwe_filename = [kwik_path filesep kwe_dir(1).name];

settings_dir = dir([kwik_path filesep '*.xml']);
settings_filename = [kwik_path filesep settings_dir(1).name];

% Get index of electrophysiology channels in recordings
kwik_settings = xml2struct(settings_filename);
ephys_ch = find(cellfun(@(x) any(strfind(x.Attributes.name,'CH')), ...
    kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL));
n_chan = length(ephys_ch);

% Get sample rate and gain
rec_info_locationInKWD = '/recordings/0';
rec_info = h5info(kwd_filename,rec_info_locationInKWD);
sample_rate_idx = cellfun(@(x) strcmp(x,'sample_rate'),{rec_info.Attributes.Name});
sample_rate = double(rec_info.Attributes(sample_rate_idx).Value);

ch_gain = cellfun(@(x) str2num(x.Attributes.gain), ...
    kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL(ephys_ch));
if length(unique(ch_gain)) == 1
    ch_gain = unique(ch_gain);
else
    error('Gains set differently for different channel: no contingency for this in param file yet')
end

% Location in the H5 file where the recordings are
locationInKWD = '/recordings/0/data';
% Load header info
info = h5info(kwd_filename,locationInKWD);

%% If local_copy selected, move data to temporary local drive for loading
if local_copy
    disp('Copying data to local drive...')
    temp_path = 'C:\temp';
    temp_filename = [temp_path filesep kwd_dir(1).name];
    if ~exist(temp_path,'dir')
        mkdir(temp_path)
    end
    copyfile(kwd_filename,temp_filename);
    disp('Done');
    
    % set data filename to the local filename
    kwd_filename = temp_filename;
end

%% Save synchronization input as .mat

disp('Getting and saving sync...');

switch sync_input
    
    case 'ttl'
        % If sync was recorded as TTL through digital input
        
        kwe_channels_loc = '/event_types/TTL/events/user_data/event_channels';
        kwe_sample_stamp_loc = '/event_types/TTL/events/time_samples';
        kwe_value_loc = '/event_types/TTL/events/user_data/eventID';
        kwe_sample_rate_loc = '/recordings/0';
        
        ttl_channels = h5read(raw_kwd_filename,kwe_channels_loc);
        ttl_samplestamp = h5read(raw_kwd_filename,kwe_sample_stamp_loc);
        ttl_values = h5read(raw_kwd_filename,kwe_value_loc);
        
        ttl_sample_rate = double(h5readatt(raw_kwd_filename,kwe_sample_rate_loc,'sample_rate'));
        
        % Save seperately all TTL events that belong to the specified sync channel
        sync = struct('timestamps',cell(size(sync_channels)), ...
            'values',cell(size(sync_channels)));
        for curr_sync = 1:length(sync_channels)
            sync_events = ttl_channels == (sync_channels(curr_sync)-1);
            sync(curr_sync).timestamps = double(ttl_samplestamp(sync_events))/ttl_sample_rate;
            sync(curr_sync).values = logical(ttl_values(sync_events));
        end
        
        sync_save_filename = [save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
    case 'adc'
        % If sync was recorded through analog input
        
        sync = struct('timestamps',cell(size(sync_channels)), ...
            'values',cell(size(sync_channels)));
        
        for curr_sync = 1:length(sync_channels)
            
            % Find the sync ADC channel index
            recorded_sync_ch = find(cellfun(@(x) strcmp(['ADC' num2str(sync_channels(curr_sync))],x.Attributes.name), ...
                kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL));
            
            % Load in the corresponding trace
            sync_trace = h5read(kwd_filename, locationInKWD,[recorded_sync_ch,1],[1,Inf]);
            
            % Binarize the sync trace based on half-max
            sync_trace = sync_trace >= max(sync_trace/2);
            
            sync_samplestamp = find((~sync_trace(1:end-1) & sync_trace(2:end)) | ...
                (sync_trace(1:end-1) & ~sync_trace(2:end)))+1;
            
            sync(curr_sync).timestamps = sync_samplestamp/sample_rate;
            sync(curr_sync).values = sync_trace(sync_samplestamp);
            
        end
        
        sync_save_filename = [save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
end

%% Filter and save spikes and LFP

% Load in each channel separately (because probably too large to load in at
% once, can be faster if multiple but easier to just to all channels
% individually for now)
% NOTE: this involves a lot of tricks in writing the file, because it
% writes in column order but should be written in row-order. This means
% that the first sample for each channel is written separately, then all
% other samples are written after n delayed bytes (and int16 = 2 bytes)

%%%% Could alternatively make downsampled LFP signal first and then load
%%%% back in and filter, but stuff from higher frequencies would be aliased
%%%% in then

disp('Loading channels, preprocessing, saving...');

if local_copy
    lfp_save_filename = [temp_path filesep 'lfp.dat'];
    spikes_save_filename = [temp_path filesep 'spikes.dat'];
else
    lfp_save_filename = [save_path filesep 'lfp.dat'];
    spikes_save_filename = [save_path filesep 'spikes.dat'];
end

lfp_fid = fopen(lfp_save_filename, 'w');
spikes_fid = fopen(spikes_save_filename, 'w');

% Figure out memory situation
dmem = memory; 
memToLeaveFree = 4 * 2^30; % num of GB to keep free
memToAllocate = dmem.MemAvailableAllArrays - memToLeaveFree;
memToAllocate = max(0, memToAllocate);
nint16s = memToAllocate/2;

n_samples = info.Dataspace.Size(2);
load_chans_n = floor(nint16s/n_samples);
%%%%% AT THE MOMENT THE FILE WRITING ISNT SET UP FOR CHUNKS: SO JUST DO ALL
load_chans_n = n_chan;
%%%%%
split_chans_n = [0,ceil(linspace(floor(n_chan/load_chans_n)-1,n_chan,floor(n_chan/load_chans_n)))];
load_chans = mat2cell(ephys_ch,1,diff(split_chans_n));

for curr_chan_chunk = 1:length(load_chans)
    
    curr_chans = load_chans{curr_chan_chunk};
    
    disp(['Channels ' num2str(curr_chans(1)) '-' num2str(curr_chans(end))]);
    
    curr_dat = h5read(kwd_filename,locationInKWD,[curr_chans(1),1],[length(curr_chans),Inf]);
    
    
    %%%%%%%%%%% TESTING
    
    if subtract_light
        
        disp('Subtracting mean-estimated light artifacts...');
        
        % get response to light onset/offset for each light
        exposure_time = mean([median(diff(sync(3).timestamps)),median(diff(sync(4).timestamps))]);
        surround_samples = -1:(floor(exposure_time*sample_rate)-1);
        
        % Define light sync channels
        light_sync = [3,4];
        
        % Go through illum 1/2 onset/offset, subtract
        for curr_sync = light_sync
            for curr_state = 0:1
                for curr_chan = 1:size(curr_dat,1)
                    curr_chan_dat = curr_dat(curr_chan,:);
                    % Define and pull samples
                    pull_samples = bsxfun(@plus,round(sample_rate* ....
                        sync(curr_sync).timestamps(sync(curr_sync).values == curr_state))',surround_samples);
                    curr_samples = curr_chan_dat(pull_samples);
                    % Get the mean light response, normalized to pre-light
                    curr_artifact = int16(nanmean(bsxfun(@minus,curr_samples,curr_samples(:,surround_samples == -1)),1));
                    curr_subtract = bsxfun(@minus,curr_samples,curr_artifact);
                    % Replace values
                    curr_chan_dat(pull_samples) = curr_subtract;
                    curr_dat(curr_chan,:) = curr_chan_dat;
                end
            end
        end
        
    end
       
    %%%%%%%%%%% TESTING
    
    
    
    % Split the signal in two:
    
    % 1) LFP (low-pass filtered and downsampled)
    disp('Filtering, downsampling, saving LFP...');
    fprintf('%2d',0);
    lfp_cutoff = 500;
    [b, a] = butter(3, lfp_cutoff/sample_rate, 'low');
    lfp_downsamp = (sample_rate/lfp_cutoff)/2;
    % filter by channel individually to avoid increasing data size
    dat_lfp = zeros(length(curr_chans),ceil(n_samples/lfp_downsamp),'int16');
    for filt_chan = 1:length(curr_chans)
        curr_lfp = filter(b,a,single(curr_dat(filt_chan,:)'))';
        curr_lfp = curr_lfp(1:lfp_downsamp:end);
        dat_lfp(filt_chan,:) = int16(curr_lfp);
        clear curr_lfp
        fprintf('%c%c%2d',8,8,floor(100*filt_chan/length(curr_chans)));
    end
    fwrite(lfp_fid,dat_lfp,'int16');
    
    %%% Use this stuff later for chunking
    %fseek(lfp_fid,(curr_chan(1)-1)*2,'bof');
    %fwrite(lfp_fid,dat_lfp(:,1),'int16');
    %fwrite(lfp_fid,dat_lfp(:,2:end),'int16',(n_chan-1)*2);
    clear dat_lfp
    
    % 2) Spikes with median across channels subtracted
    disp('Subtracting common median and saving spikes...')
    fprintf('%2d',0);
    
    % (do this in chunks: median can destroy memory)
    max_t = 1e6;
    n_chunks = ceil(size(curr_dat,2)/max_t);
    t_chunks = round(linspace(1,size(curr_dat,2),n_chunks+1));
    for curr_chunk = 1:n_chunks
        curr_t = t_chunks(curr_chunk):t_chunks(curr_chunk+1);
        % Get channel offset by median across time
        curr_chunk_offset = median(curr_dat(:,curr_t),2);
        % Remove channel offset for chunk
        curr_dat(:,curr_t) = bsxfun(@minus,curr_dat(:,curr_t),curr_chunk_offset);
        % Get and remove median across channels
        curr_chunk_median = median(curr_dat(:,curr_t),1);
        curr_dat(:,curr_t) = bsxfun(@minus,curr_dat(:,curr_t),curr_chunk_median);
        % Add back offset so there aren't jumps if median drifts
        curr_dat(:,curr_t) = bsxfun(@plus,curr_dat(:,curr_t),curr_chunk_offset);
        fprintf('%c%c%2d',8,8,floor(100*curr_chunk/n_chunks));
    end
    fwrite(spikes_fid,curr_dat,'int16');
   
    
    %%% Use this stuff later for chunking
    %fseek(spikes_fid,(curr_chan(1)-1)*2,'bof');
    %fwrite(spikes_fid,curr_dat(:,1),'int16');
    %fwrite(spikes_fid,curr_dat(:,2:end),'int16',(n_chan-1)*2);
          
end

fclose(lfp_fid);
fclose(spikes_fid);

disp('Done');

%% Save parameters/header information in separate file

params = {'raw_path',['''' kwik_path '''']; ...
    'n_channels',num2str(length(ephys_ch)); ...
    'sample_rate',num2str(sample_rate); ...
    'gain',num2str(ch_gain); ...
    'lfp_cutoff',num2str(lfp_cutoff)};

param_filename = [save_path filesep 'dat_params.txt'];

formatSpec = '%s = %s\r\n';
fid = fopen(param_filename,'w');
for curr_param = 1:size(params,1)
    fprintf(fid,formatSpec,params{curr_param,:});
end
fclose(fid);

%% If local copy, move local processed to server, delete local data copy

if local_copy    
    disp('Moving local copy to server...')
    movefile(lfp_save_filename,[save_path filesep 'lfp.dat']);
    movefile(spikes_save_filename,[save_path filesep 'spikes.dat']);
    delete(temp_filename)
    disp('Done')
end








