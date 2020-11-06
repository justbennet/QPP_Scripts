function PLTUSD2(PLc,lm,xlb,ylb)
hPL=floor((PLc(end)-PLc(1)+1)/2);
xtk=[PLc(1) PLc(1)+hPL-1 PLc(end)]; 
hold on; xticks(xtk); 
ylim(lm); ytk=lm; yticks(ytk); grid on; box on
xticklabels([]); if nargin>2 && xlb 
    a=xtk-hPL; xticklabels({num2str(a(1)),'',num2str(a(end))}); end
yticklabels([]); if nargin>3 && ylb, yticklabels(ytk); end