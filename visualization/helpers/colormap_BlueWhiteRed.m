function cm = colormap_BlueWhiteRed(n, gamma)

if nargin<1; n = 100; end
if nargin<2; gamma = 0.6; end

cm = ([n*ones(1,n), n:-1:0 ; ...
      0:n, n-1:-1:0; ...
      0:n, ones(1,n)*n]' / n).^gamma;
  
cm = cm(end:-1:1,:);  
  
colormap(cm);
