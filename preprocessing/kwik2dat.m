function kwik2dat(kwik_path,save_path,sync_channel,sync_input)
% kwik2dat(kwik_path,save_path,sync_channel,sync_input)
%
% Converts kwik format recorded from open-ephys into flat binary
% Splits signal in two: one 500 Hz low-pass LFP(lfp.dat), 
% one cross-channel median subtracted spike (spikes.dat)
% Also saves sync channel and parameters/header info
%
% kwk_path - path with kwik files
% save_path - path to save the converted data
%
% sync_channel - TTL channel with synchronization input to save
% (NOTE: input TTL channel as 1-indexed like in GUI, not 0-indexed like in
% saved data)
% sync_input - 'ttl' or 'adc', depending on digital/analog in of sync


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
ephys_ch = cellfun(@(x) any(strfind(x.Attributes.name,'CH')), ...
    kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL);

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


%% Load in the continuous electrophysiology data and save as .dat
locationInKWD = '/recordings/0/data';
info = h5info(kwd_filename,locationInKWD);
dat = h5read(kwd_filename, locationInKWD);

% Split the signal in two:

% 1) LFP (low-pass filtered and downsampled)
lfp_cutoff = 500;
[b, a] = butter(3, lfp_cutoff/sample_rate, 'low');
dat_lfp = filter(b,a,single(dat(ephys_ch,:)'))';
dat_lfp = dat_lfp(:,1:(sample_rate/lfp_cutoff)/2:end);

lfp_save_filename = [save_path filesep 'lfp.dat'];
fid = fopen(lfp_save_filename, 'w');
fwrite(fid, dat_lfp, 'int16');
fclose(fid);

% 2) Spikes with median across channels subtracted
dat_car = bsxfun(@minus,dat(ephys_ch,:),int16(median(dat(ephys_ch,:),2)));
dat_car = bsxfun(@minus,dat_car,int16(median(dat_car,1)));

spikes_save_filename = [save_path filesep 'spikes.dat'];
fid = fopen(spikes_save_filename, 'w');
fwrite(fid, dat_car, 'int16');
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
        
        sync_save_filename = [save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
    case 'adc'
        % If sync was recorded through analog input
        
        % Find the sync ADC channel index
        recorded_sync_ch = cellfun(@(x) strcmp(['ADC' num2str(sync_channel)],x.Attributes.name), ...
            kwik_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.CHANNEL_INFO.CHANNEL);
        
        sync_trace = dat(recorded_sync_ch,:) <= max(dat(recorded_sync_ch,:)/2);
        sync_samplestamp = find((~sync_trace(1:end-1) & sync_trace(2:end)) | ...
            (sync_trace(1:end-1) & ~sync_trace(2:end)));
        
        sync = struct('timestamps',[],'values',[]);
        sync.timestamps = sync_samplestamp/sample_rate;
        sync.values = sync_trace(sync_samplestamp);
        
        sync_save_filename = [save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
end

%% Save parameters/header information in separate file

params = {'raw_path',['''' kwik_path '''']; ...
    'n_channels',num2str(sum(ephys_ch)); ...
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











