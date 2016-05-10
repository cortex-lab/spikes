function kwik2dat(kwik_path,sync_channel,sync_input)
% kwik2dat(kwk_path,sync_channel)
%
% Converts kwik format recorded from open-ephys into flat binary
% also saves sync channel (currently just through TTL)
% kwk_path - path with kwik files
% sync_channel - TTL channel with synchronization input to save
% sync_input - 'ttl' or 'adc', depending on digital/analog in of sync

% Kwik files from open ephys: 
% .raw.kwd: continuous data from all channels
% .kwe: event data

%% Get filenames for each experiment in path
kwd_dir = dir([kwik_path filesep '*.kwd']);
kwd_filename = [kwik_path filesep kwd_dir(1).name];

kwe_dir = dir([kwik_path filesep '*.kwe']);
kwe_filename = [kwik_path filesep kwe_dir(1).name];

settings_dir = dir([kwik_path filesep '*.xml']);
settings_filename = [kwik_path filesep settings_dir(1).name];

% Get index of electrophysiology channels in recordings
kwik_settings = xml2struct(settings_filename);
ephys_ch = cellfun(@(x) any(strfind(x.Attributes.name,'CH')), ...
    kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL);

% Get sample rate
rec_info_locationInKWD = '/recordings/0';
rec_info = h5info(kwd_filename,rec_info_locationInKWD);
sample_rate_idx = cellfun(@(x) strcmp(x,'sample_rate'),{rec_info.Attributes.Name});
sample_rate = double(rec_info.Attributes(sample_rate_idx).Value);

%% Load in the continuous electrophysiology data and save as .dat
locationInKWD = '/recordings/0/data';
info = h5info(kwd_filename,locationInKWD);
dat = h5read(kwd_filename, locationInKWD);

ephys_save_filename = [kwik_path filesep 'ephys.dat'];
fid = fopen(ephys_save_filename, 'w');
fwrite(fid, dat(ephys_ch,:), 'int16');
fclose(fid);

%% Save synchronization input as .mat

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
        sync = struct('timestamps',[],'values',[]);
        sync_events = ttl_channels == (sync_channel-1);
        sync.timestamps = double(ttl_samplestamp(sync_events))/ttl_sample_rate;
        sync.values = logical(ttl_values(sync_events));
        
        sync_save_filename = [kwik_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
    case 'adc'
        % If sync was recorded through analog input
        
        % Find the desired ADC channel index
        recorded_sync_ch = cellfun(@(x) strcmp(['ADC' num2str(sync_channel)],x.Attributes.name), ...
            kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL);
        
        sync_trace = dat(recorded_sync_ch,:) <= max(dat(recorded_sync_ch,:)/2);
        sync_samplestamp = find((~sync_trace(1:end-1) & sync_trace(2:end)) | ...
            (sync_trace(1:end-1) & ~sync_trace(2:end)));
        
        sync = struct('timestamps',[],'values',[]);
        sync.timestamps = sync_samplestamp/sample_rate;
        sync.values = sync_trace(sync_samplestamp);
        
        sync_save_filename = [kwik_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
end





