

function d = findMinDiffs(A, B)
%function d = findMinDiffs(A, B)
% from https://uk.mathworks.com/matlabcentral/answers/79466-find-minimum-difference-between-matrices
% 
% finds the entry in B which is closest to each entry in A and return the
% distance

B = sort(B(:));     %make B vector and sort for binning convenience
edges = B(1:(end-1)) + 0.5*diff(B);   %"edges" are halfway between each B span
edges= [-inf; edges; inf];   %add +/- inf at ends to include all values
[~,whichBin] = histc(A, edges);   %get which bin contains each A
% for each value in A, the value in B which is closest to A 
%will fall between bin edges k and k+1 (i.e., k-th element of B)
%third matrix is difference between A and the B for that bin
d = A - B(whichBin);   