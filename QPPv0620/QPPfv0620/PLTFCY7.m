function PLTFCY7(c,clm,clb,xlb,ylb)
im=imagesc(c,[-1 1]); if nargin>1, im=imagesc(c,clm); end
axis square; colormap('jet');
set(im,'AlphaData',~isnan(c)); set(gca,'Color','w'); 
if nargin>2, if clb, a=colorbar; 
        ctk=[clm(1) clm(1)/2 0 clm(2)/2 clm(2)]; a.Ticks=ctk; end; end
load('myGlssr.mat','YLB','nY');
xticks(1:nY); xticklabels([]); if nargin>3 && xlb, xticklabels(YLB); end
yticks(1:nY); yticklabels([]); if nargin>4 && ylb, yticklabels(YLB); end