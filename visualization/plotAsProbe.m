

function plotAsProbe(data, xc, yc, cm, sqSizeX, sqSizeY)
% data is vector, same lenght as xc and yc
% cm is colormap, size [n x 3]
%
% other usage: if data is empty and cm's first dimension is the same length
% as xc and yc, then use those colors literally

sqCoordsX = [-sqSizeX/2 sqSizeX/2 sqSizeX/2 -sqSizeX/2];
sqCoordsY = [sqSizeY/2 sqSizeY/2 -sqSizeY/2 -sqSizeY/2];

if isempty(data)&&size(cm,1)==numel(xc)
    useLiteral = true; 
else 
    useLiteral = false;
end

if ~useLiteral
    mx = max(data); mn = min(data);
    normData = (data-mn)./(mx-mn);
    cmInds = ceil(normData*size(cm,1));
end

for q = 1:length(xc)
    
    x = xc(q)+sqCoordsX;
    y = yc(q)+sqCoordsY;        

    if useLiteral
        fill(x,y,cm(q,:), 'EdgeAlpha', 0);
    else
        fill(x,y,cmInds(q), 'EdgeAlpha', 0);
    end
    
    hold on;
end


