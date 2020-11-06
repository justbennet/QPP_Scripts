function PLTUSD(PLc,lm,xlb,ylb)
hPL=floor((PLc(end)-PLc(1)+1)/2);
xtk=[PLc(1) PLc(1)+hPL-1 PLc(end)]; 
hold on; xticks(xtk); for i=xtk, plot([i i],lm,'k'); end
ylim(lm); ytk=[lm(1) -0.2:0.1:0.2 lm(2)]; yticks(ytk); grid on; box on
plot([0 PLc(end)+hPL],[0 0],'k');
xticklabels([]); if nargin>2 && xlb 
    a=xtk-hPL; xticklabels({num2str(a(1)),'',num2str(a(end))}); end
yticklabels([]); if nargin>3 && ylb 
    yticklabels({num2str(lm(1)),'-0.2','','','','0.2',num2str(lm(2))}); end