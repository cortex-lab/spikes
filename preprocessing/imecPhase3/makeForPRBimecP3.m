

function [chanMap, chanMap0ind, xcoords, ycoords, connected, shankInd] = makeForPRBimecP3(opt)

if opt==4
    nCh = 276;
    refList = [37 76 113 152 189 228 265];
else
    nCh = 384;
    refList = [37 76 113 152 189 228 265 304 341 380];
end

chanMap = (1:nCh)';
chanMap0ind = chanMap-1;
connected = true(size(chanMap));
connected(refList) = false;
shankInd = ones(size(chanMap));
ycoords = 20*ceil(chanMap/2);
xcoords = zeros(size(chanMap));
xcoords(mod(1:nCh,4)==0) = 29;
xcoords(mod(1:nCh,4)==1) = 41;
xcoords(mod(1:nCh,4)==2) = 11;
xcoords(mod(1:nCh,4)==3) = 59;