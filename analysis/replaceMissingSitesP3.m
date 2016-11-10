

function [newData, newXC, newYC] = replaceMissingSitesP3(data, xc, yc)

sigma = 20; 

data = data(:);

missingXC = repmat([43 27]', 5, 1);
missingYC = [1:10]'*380;

missingData = zeros(size(missingXC));
for m = 1:length(missingXC)
    
    x = missingXC(m); y = missingYC(m);
    
    dists = sqrt((x-xc).^2+(y-yc).^2);
    
    % evaluate guassian of given sigma at those distances
    gaussScale = exp(-dists.^2/2/sigma^2); 
    
    gaussScale = gaussScale./sum(gaussScale); % should sum to 1 so we get the right units
    
    missingData(m) = sum(data.*gaussScale);
    
end

newData = [data; missingData];
newXC = [xc; missingXC];
newYC = [yc; missingYC];

dd = [newYC newXC newData];

dd = sortrows(dd,[1 2]);

newData = dd(:,3);
newXC = dd(:,2);
newYC = dd(:,1);