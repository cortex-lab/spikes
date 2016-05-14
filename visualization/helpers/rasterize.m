

function [xOut, yOut] = rasterize(timeStamps, varargin)
%function [xOut, yOut] = rasterize(timeStamps[, minVal, maxVal])

if ~isempty(varargin)
    minVal = varargin{1};
    maxVal = varargin{2};
else
    minVal = 0;
    maxVal = 1;
end

xOut = nan(1,length(timeStamps)*3);
xOut(1:3:end) = timeStamps; 
xOut(2:3:end) = timeStamps;

yOut = minVal + zeros(1,length(timeStamps)*3);
yOut(2:3:end) = maxVal;