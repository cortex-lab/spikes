

function splitGLXBinFile(filename, nChansTotal, subNames, subRanges, downsampleFactors, varargin)

% for opt4:
% splitGLXBinFile(filename,553,{'AP', 'LFP', 'sync'}, {1:276,277:552,553}, [1 10 10])
% for opt1,2,3:
% splitGLXBinFile(filename,769,{'AP', 'LFP', 'sync'}, {1:384,385:768,769}, [1 10 10])
%
% Exclude the '.bin' extension on filenames

if ~isempty(varargin)
    filenameOut = varargin{1};
else
    filenameOut = filename;
end

chunkSize = 1000000;

fid = []; fidOut = [];

d = dir([filename '.bin']);
nSampsTotal = d.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);
try

fid = fopen([filename '.bin'], 'r');

for n = 1:length(subNames)
    fidOut(n) = fopen([filenameOut '_' subNames{n} '.bin'], 'w');
end

theseInds = 0; chunkInd = 1;
while 1
    
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
    
    dat = fread(fid, [nChansTotal chunkSize], '*int16');
    
    if ~isempty(dat)
    
        theseInds = theseInds(end):theseInds(end)+chunkSize-1;

        for n = 1:length(subNames)
            
            fwrite(fidOut(n), dat(subRanges{n},mod(theseInds(1:size(dat,2)), downsampleFactors(n))==0), 'int16');
            
        end
    else 
        break;
    end
    
    chunkInd = chunkInd+1;
end

catch me
    
    if ~isempty(fid)
        fclose(fid);
    end
    for n = 1:length(subNames)
        if ~isempty(fidOut(n))
            fclose(fidOut(n));
        end
    end
    
    rethrow(me)
    
end
