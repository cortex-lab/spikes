

function [flipTimes, flipsUp, flipsDown] = schmittTimes(t, sig, thresh)
% function [flipTimes, flipsUp, flipsDown] = schmittTimes(t, sig, thresh)
%
% thresh is [low high]

t = t(:); % make column
sig = sig(:);

schmittSig = schmitt(sig, thresh);

flipsDown = t(schmittSig(1:end-1)==1 & schmittSig(2:end)==-1);
flipsUp = t(schmittSig(1:end-1)==-1 & schmittSig(2:end)==1);

flipTimes = sort([flipsUp; flipsDown]);