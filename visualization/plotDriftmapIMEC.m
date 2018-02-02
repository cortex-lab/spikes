% Inputs: ksDir - data directory with IMEC probe sites file forPRBimecP3opt3.mat
%         spikeTimes, spikeAmps, spikeYpos, spikeSites - names self explanatory
%         opt - optional, empty by default; 'mark' - will mark detected drifts, 'show' - will generate a different plot, 
%               where only large spikes are used, and the detection of drift locations is demonstrated
function plotDriftmapIMEC(ksDir, spikeTimes, spikeAmps, spikeYpos, spikeSites, opt)

if nargin < 6
  opt = '';
end

% LOAD SITE COORDINATES FILE
setup = load([ksDir filesep 'forPRBimecP3opt3.mat']);
sites = setup.chanMap;
connectedSites = setup.connected;
xcoords = setup.xcoords;
ycoords = setup.ycoords;

% GENERATE ELECTRODE MAP
vDim = round(numel(sites)/2);
hDim = 4;
basicMotive = [1 0 1 0 0 1 0 1];
eMapExt = repmat(basicMotive,1,vDim/2);
siteIDs = sites;
siteIDs(~connectedSites) = NaN;
eMap1(logical(eMapExt)) = siteIDs';
eMap1(~logical(eMapExt)) = NaN;
eMap1 = reshape(eMap1,hDim,vDim);
eMap1 = rot90(eMap1',2);
eMap2 = eMap1(:,[2 1 4 3]);
eMap3 = fliplr(eMap1);
eMap4 = fliplr(eMap2);
eMap = eMap2;

% PLOT DATA FOR THE WHOLE ELECTRODE
figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
plotDriftmap(spikeTimes, spikeAmps, spikeYpos, opt)
title('eMap(:,:)');

% PLOT DATA FOR EACH ELECTRODE COLUMN SEPARATELY
for c = 1:size(eMap,2)
  columnSites = eMap(:,c);
  columnSites = columnSites(~isnan(columnSites));
  spikeSitesC = zeros(size(spikeSites));
  for i = 1:numel(spikeSites)
    if any(spikeSites(i) == columnSites)
      spikeSitesC(i) = 1;
    end
  end
  spikeSitesC = logical(spikeSitesC);
  spikeTimesC = spikeTimes(spikeSitesC);
  spikeAmpsC = spikeAmps(spikeSitesC);
  spikeYposC = spikeYpos(spikeSitesC);
  figure('units', 'normalized', 'position', [0.002, .04, 1, .88]);
  plotDriftmap(spikeTimesC, spikeAmpsC, spikeYposC, opt)
  title(['eMap(:,' num2str(c) ')']);
end
