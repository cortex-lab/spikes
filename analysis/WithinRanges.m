function out = WithinRanges(x, Ranges, RangeLabel, Mode)

% out = WithinRanges(x, Ranges, RangeLabel, Mode)
% detects which points of the input vector lie within
% one of the ranges specified in the nx2 array ranges
% returns an array the size of x with a 1 if the corresponding
% point is in the ranges.
%
% ranges are (start1 stop1 ; start2 stop2 ; etc.)
% The ranges may be optionally labeled 1..nLabels
% in which case out is a matrix with one column per
% range label
%
% endpoint behaviour is inclusive
% (if i am right that sort() leaves equal in the order they
% were given in).
%
% if Mode is 'matrix' (default) it will give a matrix output
% with 1 if the point belongs to that range label
% if Mode is 'vector' it will give a vector, specifying the range
% of each point (gives error if any point belongs to more than 1)

% reshape x to a vector
x = x(:);

% get size info
nPoints = length(x);
nRanges = size(Ranges,1);

if nargin<3
	RangeLabel = ones(nRanges, 1);
end
nLabels = max(RangeLabel);
if nargin<4
    Mode = 'matrix';
end


if nRanges==0
    out=zeros(size(x));
    return;
end

% check End comes after Start in each case
if any(Ranges(:,2)<Ranges(:,1))
    error('End should come after Start!');
end

% make array containing points, starts and finishes

%ToSort = [x ; Ranges(:,1) ; Ranges(:,2)];
ToSort = [Ranges(:,1) ; x ; Ranges(:,2)]; % this order means it will be inclusive
% sort it
[Sorted Index] = sort(ToSort);

% Make delta array containing 1 for every start and -1 for every stop
% with one column for each range label

if strcmp(Mode, 'matrix')
    Delta = zeros(nPoints+2*nRanges,nLabels);
    Delta(sub2ind([nPoints+2*nRanges,nLabels],1:nRanges,RangeLabel')) = 1;
    Delta(sub2ind([nPoints+2*nRanges,nLabels],nPoints+nRanges+(1:nRanges),RangeLabel')) = -1;

    %Arrange it in order
    DeltaSorted = Delta(Index,:);

    % take cumulative sums
    Summed = cumsum(DeltaSorted);

	% and reorder back to the original order
	ReOrdered = zeros(nPoints+2*nRanges,nLabels);
	ReOrdered(Index,:) = Summed;
	
	out = ReOrdered(nRanges+1:nPoints+nRanges,:);
elseif strcmp(Mode, 'vector')
    nDelta = zeros(nPoints+2*nRanges,1);
    nDelta(1:nRanges) = 1;
    nDelta(nPoints+nRanges+(1:nRanges)) = -1;
    rDelta = zeros(nPoints+2*nRanges,1);
    rDelta(1:nRanges) = RangeLabel;
    rDelta(nPoints+nRanges+(1:nRanges)) = -RangeLabel;
    
    nDeltaSorted = nDelta(Index);
    rDeltaSorted = rDelta(Index);

    % take cumulative sums
    nSummed = cumsum(nDeltaSorted);
    rSummed = cumsum(rDeltaSorted);

	% and reorder back to the original order
	nReOrdered(Index) = nSummed;
	rReOrdered(Index) = rSummed;
    
    if any(nReOrdered(nRanges+1:nPoints+nRanges)>1)
        error('Some points belong to more than one range');
    else
        out = rReOrdered(nRanges+1:nPoints+nRanges);
    end
end
    




return
