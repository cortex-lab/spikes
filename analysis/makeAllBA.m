

function [allba, bins] = makeAllBA(st, clu, eventTimes, binSize, win)

binBorders = win(1):binSize:win(2);
numBins = length(binBorders)-1;
bins = binBorders(1:end-1)+binSize/2;

% now using histdiffMulti. 
clu = clu-min(clu)+1; % make sure they are all positive
cIDs = unique(clu); nC = numel(cIDs);

clu = clu(st>min(eventTimes)+win(1) & st<max(eventTimes)+win(2));
st = st(st>min(eventTimes)+win(1) & st<max(eventTimes)+win(2)); % cut out all distant st to start

% first thing is to replace clu with a vector that has the same size but
% contains (the indices of cids)-1 
inds = nan(max(cIDs),1);
inds(cIDs) = 0:nC-1; % these are the indices we want them to have
g = inds(clu);

[st,ii] = sort(st);
g = g(ii);


% replace NaN event times with interpolated values so the array is all
% non-nan and still sorted. then we delete those rows later.
nanevts = isnan(eventTimes); 
eventTimes(nanevts) = interp1(find(~nanevts), eventTimes(~nanevts), find(nanevts),'linear', 'extrap');

if isempty(st)
    n = zeros(1, prod([numel(bins) numel(eventTimes) nC]));
else
    [n,~] = histdiffMulti(st, eventTimes, binBorders, g, nC);
end

outputSize = [numel(bins) numel(eventTimes) nC];
if numel(n)~=prod(outputSize)
    fprintf(1, 'histdiffmulti: wrong size returned, got %d elements, wanted [%d %d %d]=%d\n',...
        numel(n), outputSize(1), outputSize(2), outputSize(3), prod(outputSize));
end
n = reshape(n, outputSize);
allba = permute(n, [3 2 1]); % returned order is cluster x evNum x timebins

% remove the nan events
allba(:,nanevts,:) = NaN;

% ** OLD METHOD, v. slow
% cIDs = unique(clu);
% allba = zeros(numel(cIDs), numel(eventTimes), numBins);
% 
% clu = clu(st>min(eventTimes)+win(1) & st<max(eventTimes)+win(2));
% st = st(st>min(eventTimes)+win(1) & st<max(eventTimes)+win(2)); % cut out all distant st to start
% 
% % there will be a faster way to implement this with WithinRanges, at least
% % in the case of non-overlapping ranges
% % theseST = st(clu==cIDs(c));
% theseST = arrayfun(@(x)st(clu==x), cIDs, 'uni', false);
% for c = 1:numel(cIDs)
% %     fprintf(1, '%d/%d\n', c, numel(cIDs))
%     
% 
%     if ~isempty(theseST{c})
%         [ba, bins] = timestampsToBinned(theseST{c}, eventTimes, binSize, win);
% 
%         allba(c, :,:) = ba;
%     end
% end
% 
% bins = binBorders(1:end-1)+binSize/2; 
