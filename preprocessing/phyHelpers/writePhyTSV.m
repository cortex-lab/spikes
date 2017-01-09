

function writePhyTSV(ksDir, labelName, cids, vals)
% function writePhyTSV(ksDir, labelName, cids, vals)
% Writes a .tsv file that phy can read to apply a new sortable column in cluster
% view. E.g. calculate the amplitude of every template and use this to
% write a .tsv file that will allow you to sort by that amplitude. 
%
% - ksDir - the directory of your kilosort results (i.e. phy files)
% - labelName - a string, what to call the label (e.g. 'amplitude')
% - cids - the cluster ID numbers, a vector, all integers
% - vals - can be a cell array of strings, or a vector of numbers. Should be
% the same size as cids. 

fid = fopen(fullfile(ksDir, ['cluster_' labelName '.tsv']), 'w');

fwrite(fid, sprintf('cluster_id\t%s\r\n', labelName));

[cids, ii] = sort(cids);
vals = vals(ii);

for c = 1:length(cids)
    if iscell(vals)
        fwrite(fid, sprintf('%d\t%s\r\n', cids(c), vals{c}));
    else
        fwrite(fid, sprintf('%d\t%s\r\n', cids(c), num2str(vals(c))));
    end
end
fclose(fid);