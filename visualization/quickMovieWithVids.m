

function quickMovieSpikes(mouseName, thisDate, tag, type)

% determine which timeline to use


% load spikes


expPath = dat.expPath(mouseName, thisDate, expNum, 'main', 'master');

% load traces
load(dat.expFilePath(mouseName, thisDate, expNum, 'Timeline', 'master'));
traces = prepareTimelineTraces(Timeline);

% load videos
auxVid = prepareAuxVids(mouseName, thisDate, expNum);

writeMovieLocation = fullfile(expPath, sprintf('spikes_%s_%s_%d', mouseName, thisDate, expNum));


switch type
    case 'mua'
        
    case 'clu'
        
end
movieWithTracesSpikes(spikeTimes, spikeClu, traces, auxVid, pars)