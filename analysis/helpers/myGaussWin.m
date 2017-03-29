
function gw = myGaussWin(stdev, Fs)
% function gw = myGaussWin(stdev, Fs)
% A gaussian window with specified stdev in units relative to the sampling
% frequency, and normalized in amplitude
stdevSamps = round(stdev*Fs);
gw = gausswin(stdevSamps*6,3);
gw = gw./sum(gw);

