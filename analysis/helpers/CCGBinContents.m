function SpikePairs = CCGBinContents(T1, T2, BinStart, BinEnd)
% SpikePairs = CCGBinContents(T1, T2, BinStart, BinEnd)
%
% tells you which pairs of spikes contributed to a CCG bin
%
% T1 and T2 give the spike times of the 2 cells.
%
% BinStart and BinEnd describe the bin in the same units.
%
% output SpikePairs is an nx2 matrix giving pairs of spikes
% in the CCG bin.

SpikePairs = zeros(0,2);

Size1 = size(T1, 1);
Size2 = size(T2, 1);

RangeStart = 1;
RangeEnd = 1;

for Sp1 = 1:Size1

  t1 = T1(Sp1);

  Sp2 = find(T2 >= t1 + BinStart & T2 <= t1+BinEnd);

  SpikePairs = [SpikePairs ; repmat(Sp1, length(Sp2), 1), Sp2(:) ];
end;

