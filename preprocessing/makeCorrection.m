
function [corrFun, b] = makeCorrection(correctTo, correctFrom, makePlots)
% function [corrFun, b] = makeCorrection(correctTo, correctFrom, makePlots)
%
% make a function that linearly predicts correctTo from correctFrom, both vectors of same
% length
% - now can attempt to fit for cases with vectors of different lengths,
% i.e. when one 
s1 = correctTo(:);
s2 = correctFrom(:);

if length(s1)==length(s2)
    b = regress(s1, [s2 ones(size(s2))]);
else
    fprintf(1, 'warning: may be missing events in one of the two inputs! trying to fit anyway\n')
    op = optimoptions('fmincon', 'TolX', 1e-9);
    %b = fminunc(@(xIn)totalMinDiffs(xIn, s1, s2), [1. s1(1)-s2(1)]', op);
    
    % using the constrained version allows us to exclude the trivial
    % solution of b(1) = 0, i.e. multiplying every entry by zero gets them
    % all very close to one of the others, trivially satisfying the
    % minimization of distances.
    b = fmincon(@(xIn)totalMinDiffs(xIn, s1, s2), [1. s1(1)-s2(1)]', [], [], [], [], [0.5 -inf], [1.5 inf], [], op);
end

corrFun = @(s2)[s2(:) ones(size(s2(:)))]*b;

if makePlots
    figure; 

    startOffset = s1(1)-s2(1);
    
    plot(s1, findMinDiffs(s1,corrFun(s2)), '.');


    hold on; 
    plot(corrFun(s2), corrFun(s2)-s2-startOffset, 'r');

    xlabel('time in To system(sec)');
    ylabel('difference between correctTo and correctFrom (sec)')
    legend({'events measured on each system', 'best fit'});
    % makepretty;
end

fprintf(1, 'drift rate was %.3f msec / hour\n', (b(1)-1)*60*60*1000);


