function PLTT(T,PLc,lm,clbr,xlb,ylb)
load('myGlssr.mat','nY','ibY','iG2Y','ihY','YLB');
To=T; for i=1:nY, To(ibY(i)+1:ibY(i+1),:)=T(iG2Y{i},:); end
if nargin>2, im=imagesc(To,lm); else, im=imagesc(To,[-1 1]); end
colormap('jet'); set(im,'AlphaData',~isnan(To)); hold on; 
xtk=[PLc(1) PLc(end)-PLc(1)+1 PLc(end)];
[nX,PLe]=size(T);
for i=xtk, plot([i i],[0 nX+1],'k'); end; xticks(xtk);
for i=ibY(2:end-1), plot([0 PLe+1],[i i],'k'); end; yticks(ihY);
if nargin>3, if clbr, a=colorbar; 
        ctk=[lm(1) lm(1)/2 0 lm(2)/2 lm(2)]; a.Ticks=ctk; end; end
xticklabels([]); if nargin>4 && xlb, xticklabels(xtk-PLc(1)+1); end
yticklabels([]); if nargin>5 && ylb, yticklabels(YLB); end