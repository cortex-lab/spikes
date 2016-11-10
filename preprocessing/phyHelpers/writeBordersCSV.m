

function writeBordersCSV(borders, fn)

fid = fopen(fn, 'w');

fprintf(fid, 'area\tlowerBorder\tupperBorder\r\n');

fn = fields(borders);

for f = 1:length(fn)
    fprintf(fid, '%s\t%d\t%d\r\n', fn{f}, borders.(fn{f})(1), borders.(fn{f})(2));
end

fclose(fid);