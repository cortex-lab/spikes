
function sp = loadAllKsDir(mouseName, thisDate)

rootE = dat.expPath(mouseName, thisDate, 1, 'main', 'master');
root = fileparts(rootE);
% ksRoot = fullfile('\\basket.cortexlab.net\data\nick', mouseName, thisDate);
alignDir = fullfile(root, 'alignments');

% d = dir(fullfile(root, 'ephys*'));
% if numel(d)>1
%     for q = 1:numel(d)
%         tags{q} = d(q).name(7:end);
%     end
% end
[tags, hasEphys] = getEphysTags(mouseName, thisDate);

for q = 1:length(tags)
    fprintf(1, 'loading %s\n', tags{q});
    
    tag = tags{q};
    
%     ksDir = fullfile(ksRoot, ['ephys_' tag]);
    ksDir = getKSdir(mouseName, thisDate, tag);
    rawDir = fullfile(root, ['ephys_' tag]);
    addGainToParamsPy(ksDir, rawDir)
    
    sptemp = loadKSdir(ksDir);
    sptemp.name = tag;
    if q==1
        sp = sptemp;
    else
        
        b = readNPY(fullfile(alignDir, ...
            sprintf('correct_ephys_%s_to_ephys_%s.npy', tag, tags{1})));
        
        sptemp.st = applyCorrection(sptemp.st, b);
        
        sp(q) = sptemp;
    end
end
    
for q = 1:length(sp)
    fprintf(1, 'computing stuff for %s\n', tags{q});
    
    [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
        templatePositionsAmplitudes(sp(q).temps, sp(q).winv, sp(q).ycoords, sp(q).spikeTemplates, sp(q).tempScalingAmps);

    
    
    if isfield(sp(q), 'gain')
        spikeAmps = spikeAmps*sp(q).gain;
        tempAmps = tempAmps*sp(q).gain;
    end

    sp(q).spikeAmps = spikeAmps;
    sp(q).spikeDepths = spikeDepths;
    sp(q).templateYpos = templateYpos;
    sp(q).tempAmps = tempAmps;
    sp(q).tempsUnW = tempsUnW;
    sp(q).tempDur = tempDur;
    sp(q).tempPeakWF = tempPeakWF;
    
    
end