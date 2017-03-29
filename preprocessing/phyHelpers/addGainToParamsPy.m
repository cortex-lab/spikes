
function addGainToParamsPy(ksDir, rawDir)

paramsFn = fullfile(ksDir, 'params.py');
p = loadParamsPy(paramsFn);
if ~isfield(p, 'gain')
    [~, rawFn] = fileparts(p.dat_path);
    if strcmp(rawFn(end-3:end), '_CAR')
        metaFn = [rawFn(1:end-4) '.meta'];
    else
        metaFn = [rawFn '.meta'];
    end
    metaPath = fullfile(rawDir, metaFn);
    if exist(metaPath, 'file')
        meta = readSpikeGLXmeta(metaPath);
        if isfield(meta, 'uV_per_bit')
            writeToParamsPy(paramsFn, 'gain', meta.uV_per_bit);
        else
            fprintf(1, 'meta did not have field uV_per_bit');
        end
        if isfield(meta, 'uV_per_bit_lfp')
            writeToParamsPy(paramsFn, 'gainLFP', meta.uV_per_bit_lfp);
        end
    else
        fprintf(1, 'did not find meta file at %s\n', metaPath);
    end
else
%     fprintf(1, 'params.py at %s already has gain field\n', paramsFn);
end
