

function outData = applyCorrection(inData, b)
% function outData = applyCorrection(inData, b)
%
% align a timeseries according to correction factors ("b") made with
% companion function "makeCorrection"

outData = [inData(:) ones(size(inData(:)))]*b;
