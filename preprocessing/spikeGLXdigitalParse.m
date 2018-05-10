
function eventTimes = spikeGLXdigitalParse(digitalChannel, Fs)
% function eventTimes = spikeGLXdigitalParse(digitalChannel, Fs)
%
% returns the event times of all 16 digital inputs recorded by SpikeGLX
% They are returned in a length=16 cell array. Each cell has three cells,
% which contain the times of all events (sorted), the times of onsets, and
% the times of offsets, respectively. 
%
% Fs is optional input argument - if supplied, output will be in units of
% seconds, otherwise samples

digitalChannel = digitalChannel(:); % column

if isa(digitalChannel, 'int16')
    % change to uint16 but don't change the underlying bits
    digitalChannel = typecast(digitalChannel, 'uint16');
end

if exist('de2bi')==0
    bitRepresentation = arrayfun(@(x)str2double(x), dec2bin(digitalChannel,16));
else
    % faster but this is in communication systems toolbox
    bitRepresentation = de2bi(digitalChannel,16); % size is nSamples x 16
end

for b = 1:16
    db = diff([0;double(bitRepresentation(:,b))]);
    bOn = find(db>0); 
    bOff = find(db<0);
    
    if nargin>1
        eventTimes{b}{1} = sort([bOn; bOff])/Fs;
        eventTimes{b}{2} = bOn/Fs;
        eventTimes{b}{3} = bOff/Fs;
    else
        eventTimes{b}{1} = sort([bOn; bOff]);
        eventTimes{b}{2} = bOn;
        eventTimes{b}{3} = bOff;
    end
end


