

function quickMovieSpikesALF(mouseName, thisDate, expNum, tag, type)


root = fileparts(dat.expPath(mouseName, thisDate, 1, 'main', 'master'));
alfDir = fullfile(root, 'alf');

% load traces
traceNames = {'eye.area', 'lickSignal.trace', 'wheel.velocity'};

for tInd = 1:length(traceNames)
    tr = readNPY(fullfile(alfDir, [traceNames{tInd} '.npy']));
    base = traceNames{tInd}(1:find(traceNames{tInd}=='.')-1);
    t = readNPY(fullfile(alfDir, [base '.timestamps.npy']));
    tVec = interp1(t(:,1), t(:,2), 0:numel(tr)-1);
    traces(tInd).t = tVec;
    traces(tInd).v = tr;
    traces(tInd).name = traceNames{tInd};
    if strcmp(traceNames{tInd}, 'eye.area') 
        % normalize to just the time during the task
        trials = readNPY(fullfile(alfDir,'cwTrials.intervals.npy'));
        eyeInTr = tr(tVec>min(trials(:)) & tVec<max(trials(:)));        
        traces(tInd).lims = prctile(eyeInTr, [5 95]);
    else
        traces(tInd).lims = [-1 1]*max(abs(tr))*0.75;
    end
end

t = readNPY(fullfile(alfDir, 'cwStimOn.times.npy'));
[xx,yy] = rasterize(t);
traces(end+1).t = xx;
traces(end).v = yy;
traces(end).name = 'stimOn';

t = readNPY(fullfile(alfDir, 'cwGoCue.times.npy'));
[xx,yy] = rasterize(t);
traces(end+1).t = xx;
traces(end).v = yy;
traces(end).name = 'goCue';

t = readNPY(fullfile(alfDir, 'cwFeedback.times.npy'));
[xx,yy] = rasterize(t);
traces(end+1).t = xx;
traces(end).v = yy;
traces(end).name = 'feedback';

% load videos
auxVid = prepareAuxVids(mouseName, thisDate, expNum);
faceT = readNPY(fullfile(alfDir, 'face.timestamps.npy'));
tVec = interp1(faceT(:,1), faceT(:,2), 0:numel(auxVid(1).data{2})-1);
auxVid(1).data{2} = tVec;
eyeT = readNPY(fullfile(alfDir, 'eye.timestamps.npy'));
tVec = interp1(eyeT(:,1), eyeT(:,2), 0:numel(auxVid(2).data{2})-1);
auxVid(2).data{2} = tVec;

st = readNPY(fullfile(alfDir, tag, 'spikes.times.npy'));
depths = readNPY(fullfile(alfDir, tag, 'spikes.depths.npy'));

borders = readNPY(fullfile(alfDir, tag, 'spikes.times.npy'));

anatData.coords = [sp(spInd).xcoords sp(spInd).ycoords];
anatData.borders = borders;


% create a simulated set of "waveforms" that will just highlight the
% correct segment of the probe
uClu = unique(clu); 
fakeWF = zeros(numel(uClu), numel(sp(spInd).xcoords));
ycBins = ceil(sp(spInd).ycoords/depthBin)*depthBin;
for c = 1:numel(uClu)
    fakeWF(c,ycBins==uClu(c)) = 1;
end
anatData.wfLoc = fakeWF;

switch type
    case 'mua'
        depthBinSize = 80;
        clu = ceil(depths/depthBinSize);
        pars.smoothSizeT = 0;
        
        % create a simulated set of "waveforms" that will just highlight the
        % correct segment of the probe
        uClu = unique(clu);
        fakeWF = zeros(numel(uClu), numel(xcoords));
        ycBins = ceil(ycoords/depthBinSize)*depthBinSize;
        for c = 1:numel(uClu)
            fakeWF(c,ycBins==uClu(c)) = 1;
        end
        anatData.wfLoc = fakeWF;
    case 'clu'
        clu = readNPY(fullfile(alfDir, tag, 'spikes.clusters.npy'));
        % now re-sort the clu numbers by depth
end
pars.winSize = [-2 1];
pars.normalize = false;
pars.binSize = 0.002;
movieWithTracesSpikes(st, clu, traces, auxVid, pars)