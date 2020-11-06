function PLTFC(c,clm,clb,xlb,ylb)
im=imagesc(c,[-1 1]); if nargin>1, im=imagesc(c,clm); end
axis square; colormap('jet'); set(im,'AlphaData',~isnan(c)); 
if nargin>2 && clb, colorbar; end
load('myGlssr.mat','ibY','ihY','YLB','nG');
hold on; for i=ibY(2:end-1), plot([0 nG+1],[i i],'k'); end
for i=ibY(2:end-1), plot([i i],[0 nG+1],'k'); end % ,'Linewidth',2
xticks(ihY); xticklabels([]); if nargin>3 && xlb, xticklabels(YLB); end
yticks(ihY); yticklabels([]); if nargin>4 && ylb, yticklabels(YLB); end