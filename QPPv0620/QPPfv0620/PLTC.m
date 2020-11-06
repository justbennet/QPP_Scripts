function PLTC(c,clm,clb,xlb,ylb)
im=imagesc(c,[-1 1]); if nargin>1, im=imagesc(c,clm); end
axis image; colormap('jet');
set(im,'AlphaData',~isnan(c));
if nargin>2, if clb, colorbar; end; end
[ny,nx]=size(c);
xticks(1:nx); if nargin>3, xticklabels(xlb); end
yticks(1:ny); if nargin>4, yticklabels(ylb); end