
function S = readSpikeGLXmeta(vcFname)
% function S = readSpikeGLXmeta(vcFname)
%
% returns a struct containing the contents of a SpikeGLX-generated .meta
% file. 
%
% From James Jun, 2016-12-31

S = [];
viRef_imec3 = [37 76 113 152 189 228 265 304 341 380];

% Read file
if nargin < 1
    [FileName,PathName,FilterIndex] = uigetfile();
    vcFname = fullfile(PathName, FileName);
    if ~FilterIndex
        return; 
    end
end

try
    %Read Meta
    S = text2struct(vcFname);    
    S.vcDataType = 'int16'; %whisper standard

    %convert new fields to old fields    
    if isfield(S, 'niSampRate')        
        % SpikeGLX
        S.nChans = S.nSavedChans;
        S.sRateHz = S.niSampRate;
        S.rangeMax = S.niAiRangeMax;
        S.rangeMin = S.niAiRangeMin;
        S.auxGain = S.niMNGain;
        try
            S.outputFile = S.fileName;
            S.sha1 = S.fileSHA1;      
            S.vcProbe = 'imec2';
        catch
            S.outputFile = '';
            S.sha1 = [];      
            S.vcProbe = '';
        end
        S.ADC_bits = 16;
    elseif isfield(S, 'imSampRate')
        % IMECIII
        S.nChans = S.nSavedChans;
        S.sRateHz = S.imSampRate;
        S.rangeMax = S.imAiRangeMax;
        S.rangeMin = S.imAiRangeMin;
        S.ADC_bits = 10;  %10 bit adc but 16 bit saved
        
        
        % new, for p3b compatibility, in which the format has changed. 
        imroEntries = regexp(S.imroTbl, '[\(].*?[\)]', 'match');
        summaryData = str2num(imroEntries{1}(2:end-1));
        firstChannel = str2num(imroEntries{2}(2:end-1));
        S.auxGain = firstChannel(4);
        S.auxGain_lfp = firstChannel(5);
        S.nSites = summaryData(end);
        
%         vnIMRO = textscan(S.imroTbl, '%d', 'Delimiter', '( ),');
%         vnIMRO = vnIMRO{1};
%         S.auxGain = double(vnIMRO(9)); %hard code for now;
%         S.auxGain_lfp = double(vnIMRO(10)); %hard code for now;
%         S.vcProbe = sprintf('imec3_opt%d', vnIMRO(3));
%         S.nSites = vnIMRO(4);
%         S.viSites = setdiff(1:S.nSites, viRef_imec3); %sites saved
        try
            S.S_imec3 = imec3_imroTbl(S);
        catch
            S.S_imec3 = [];
        end
    else
        S.vcProbe = 'generic';
        S.ADC_bits = 16;
    end
    
     %number of bits of ADC [was 16 in Chongxi original]
    S.scale = ((S.rangeMax-S.rangeMin)/(2^S.ADC_bits))/S.auxGain * 1e6;  %uVolts
    S.uV_per_bit = S.scale;
    if isfield(S, 'auxGain_lfp')
        S.uV_per_bit_lfp = S.scale * S.auxGain / S.auxGain_lfp;
    end
catch
    disp(lasterr);
end
end 