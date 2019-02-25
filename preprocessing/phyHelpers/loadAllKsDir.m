
function sp = loadAllKsDir(mouseName, thisDate)
% function sp = loadAllKsDir(mouseName, thisDate)
% load all the spikes from a day into one struct, aligning them to each
% other



root = getRootDir(mouseName, thisDate);

alignDir = fullfile(root, 'alignments');

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