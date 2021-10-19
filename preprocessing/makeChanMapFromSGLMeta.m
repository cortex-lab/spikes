

function cm = makeChanMapFromSGLMeta(m, varargin)
% function cm = makeChanMapFromSGLMeta(m[, shankSep, rowSep, colSep])

if nargin>1
    shankSep = varargin{1}; 
    rowSep = varargin{2}; 
    colSep = varargin{3};
else
    % assume NP 2.0 parameters
    shankSep = 250; 
    rowSep = 15; 
    colSep = 32;
end

shankMap = m.snsShankMap; 

% why tf does this keep changing...?
if isfield(m, 'acqApLfSy')
    nCh = str2num(m.acqApLfSy);
elseif isfield(m, 'snsApLfSy')
    nCh = str2num(m.snsApLfSy);
elseif isfield(m, 'nChans')
    nCh = m.nChans-1;
elseif isfield(m, 'nSavedChans')
    nCh = m.nSavedChans-1;
end


chanMap = [1:nCh(1)]'; 
chanMap0ind = chanMap-1;
connected = true(size(chanMap)); 

openParen = find(shankMap=='('); 
closeParen = find(shankMap==')'); 
for c = 1:nCh(1)
    thisMap = shankMap(openParen(c+1)+1: closeParen(c+1)-1); 
    thisMap(thisMap==':') = ',';
    n = str2num(thisMap); 
    xcoords(c) = (n(1)-1)*shankSep + (n(2)-1)*colSep; 
    ycoords(c) = (n(3)-1)*rowSep; 
end

cm = struct();
cm.chanMap = chanMap; 
cm.chanMap0ind = chanMap0ind;
cm.xcoords = xcoords'; 
cm.ycoords = ycoords'; 
cm.connected = connected;
[~,name] = fileparts(m.imRoFile); 
cm.name = name;