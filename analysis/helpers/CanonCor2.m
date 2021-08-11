% [a b R2 V] = CanonCor2(Y, X)
%
% does a sort of canonical correlation analysis on two sets of data X and Y
% except that now it finds the linear combinations of X that predict
% the largest variance fractions of Y.
%
% You should think of Y as the dependent variable, and X as the independent
% variable.
%
% R2 is the fraction of the total variance of Y explained by
% the nth projection
%
% the approximation of Y based on the first n projections is:
% Y = X * b(:,1:n) *a'(:,1:n);
%
% the nth variable for the ith case gives the approximation
% Y(i,:)' = a(:,n) * b(:,n)' * X(i,:)'
%
%
% V is the actual value of the nth linear combination of X.

function [a, b, R2, V] = CanonCor2(Y, X)

% Make covariance matrices
XSize = size(X, 2);
%YSize = size(Y, 2);

BigCov = cov([X, Y]);
CXX = BigCov(1:XSize, 1:XSize);
%CYY = BigCov(XSize+1:end, XSize+1:end);
%CXY = BigCov(1:XSize, XSize+1:end);
CYX = BigCov(XSize+1:end, 1:XSize);

eps = 1e-7; 
CXX = CXX+eps*eye(XSize); % prevents imaginary results in some cases

CXXMH = CXX ^ -0.5; 

% matrix to do svd ...

M = CYX * CXXMH;

% do svd
[d, s, c] = svd(M, 0);

b = CXXMH * c;
a = d * s;

R2 = (diag(s).^2)/sum(var(Y));

if (nargout > 3)
	V = X*b;
end;
