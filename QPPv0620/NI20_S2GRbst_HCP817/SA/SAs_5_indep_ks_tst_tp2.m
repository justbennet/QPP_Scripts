
%% Testing significance of time differences btwn region using independent 
%% t-test and ks-test, without matching group sizes
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=p2p.name;
load([p1 p2p],'p2V','nP','nSD','nsd','PL','nvx','nrgn','irgn','nmrgn1','d2SAplt'); 
nth=15; benf=nrgn*(nrgn-1)/2; alph=0.01/benf;
load([p1 p2V],'QPPv_te','QPPv_tm'); 
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close;

%% Common part to t-/ks-test
QPPv_tp2=cell(nrgn,nP,nSD); QPPv_tp2(:)={single([])};
for ip=1:nP
for isd=1:nsd(ip)
    tp=QPPv_te(:,1,ip,isd); tm=QPPv_tm(ip,isd);
    for ig=1:nrgn
        tp2=nan(nvx,1,'single'); tp2(irgn{ig})=tp(irgn{ig});
        tp2(isnan(tp2))=[]; tp2(tp2<tm)=[]; n=length(tp2); 
        if n>=nth
            QPPv_tp2{ig,ip,isd}=tp2; 
        end
    end
end
end; clear tp tp2 tm

%% uncomment t-/ks-test whichever is to test & plot
pttst_tp2=nan(nrgn,nrgn,nP,nSD,'single');
for ip=1:nP
for isd=1:nSD 
for ig1=1:nrgn
    tp1=QPPv_tp2{ig1,ip,isd}; l1=length(tp1);
    if l1
    for ig2=1:nrgn
        tp2=QPPv_tp2{ig2,ip,isd}; l2=length(tp2);
        if l2
            [h,p]=ttest2(tp1,tp2,'Alpha',alph,'vartype','unequal');
%             [h,p]=kstest2(tp1,tp2,'Alpha',alph);
            if ~h, pttst_tp2(ig1,ig2,ip,isd)=1; 
            else, pttst_tp2(ig1,ig2,ip,isd)=p; end
        end
    end
    end
end
end
end

s=ss2; s(1)=1; s(3)=ss(3); 
for isd=1:nSD
    figure; set(gcf,'Position',s);
for ip=1:nP
    p=pttst_tp2(:,:,ip,isd);
    subplot(1,nP,ip),im=imagesc(p); axis square; set(im,'AlphaData',~isnan(p));
    xticks(1:nrgn); yticks(1:nrgn); xticklabels(nmrgn1); yticklabels(nmrgn1);
end
    saveas(gcf,[d2SAplt 'SAs_5_indep-ttst-tp2_SD' num2str(isd) '.png']); close
%     saveas(gcf,[d2SAplt 'SAs_5_ks-tst-tp2_SD' num2str(isd) '.png']); close
end