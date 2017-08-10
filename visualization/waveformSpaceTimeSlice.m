

function wfSlice = waveformSpaceTimeSlice(t, coords, slicePos)
% function wfSlice = waveformSpaceTimeSlice(t, coords, slicePos)
% 
% Interpolates a "slice" of the spatiotemporal waveform, in the x-dimension
% (so you get a Y by Time slice)

%%
outlierDistance = 50; 
spaceRange = 150; %µm, how much of the space to return - peak location +/- this
spacing = 4; %µm, resolution of the interpolation

[~,maxChan] = max(max(abs(t),[],2));
maxYpos = coords(maxChan,2);
inclChans = coords(:,2)>maxYpos-spaceRange & coords(:,2)<maxYpos+spaceRange;

minX = min(coords(:,1)); maxX = max(coords(:,1));
leftEdgeY = coords(inclChans&coords(:,1)==minX,2);
rightEdgeY = coords(inclChans&coords(:,1)==maxX,2);

nLeft = numel(leftEdgeY);
outliersLeft = [repmat(minX-outlierDistance, nLeft, 1) leftEdgeY];
nRight = numel(rightEdgeY);
outliersRight = [repmat(maxX+outlierDistance, nRight, 1) rightEdgeY];

x = [coords(inclChans,1); outliersLeft(:,1); outliersRight(:,1)];
y = [coords(inclChans,2); outliersLeft(:,2); outliersRight(:,2)];

% [icx, icy] = meshgrid(min(x):spacing:max(x), min(y):spacing:max(y));
icy = min(y):spacing:max(y);
icx = slicePos*ones(size(icy));

% make interpolant
F = scatteredInterpolant(x, y, zeros(size(x)), 'natural');

% compute the slice
nSkip = 0; %time points to skip at the beginning of the waveform, it's a kilosort thing
wfSlice = zeros(numel(icy), size(t,2)-nSkip);

%%

% I think you could do the whole interpolation at once in 3D with
% scatteredInterpolant, however I don't really want the different time
% slices interacting with each other, and I'm not sure how the different
% scaling of the time and space axes would be handled, so...
for tInd = nSkip+1:size(t,2) 
    wfData = t(inclChans,tInd)*100; % *100? can't recall why
    
    v = double([wfData; zeros(size(outliersLeft(:,1))); zeros(size(outliersRight(:,1)))]);
    
    F.Values = v;
    vq = F(icx, icy);
    wfSlice(:,tInd-nSkip) = vq;
end