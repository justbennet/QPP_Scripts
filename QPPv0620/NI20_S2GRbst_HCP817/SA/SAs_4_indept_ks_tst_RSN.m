
%% Testing significance of time differences btwn cortical RSNs using
%% independent t-test and ks-test, without matching group sizes
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=p2p.name;
load([p1 p2p],'p2V','p2u','nY','ivx','nP','nSD','nsd','PLe','PL','PLc','nvx','d2SAplt'); 
benf=nY*(nY-1)/2; alph=0.01/benf;
load([p1 p2V],'QPPv','QPPv_io','QPPv_te'); 
load([p1 p2u 'myGlssr.mat'],'Y','YLB'); Y=Y(ivx);
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close; nmte='pz'; 

%% Common part to t-/ks-test
QPPv2y_te=cell(nP,nSD); QPPv2y_te(:)={cell(nY,1)};
for ip=1:nP
for isd=1:nsd(ip)
    
    T=QPPv{ip,isd}; io=QPPv_io{ip,isd};
    ix=1:nvx; ix(io)=[]; T(ix,:)=nan;
    
    tz=nan(nvx,1,'single');
    for iv=io'
        z=T(iv,PLc); t=find(z.*circshift(z,-1)<=0); t(t==PL)=[];
        if length(t)==1, tz(iv)=t; end
    end
    
    for iy=1:nY
        tez=[QPPv_te(Y==iy,:,ip,isd) tz(Y==iy)]; 
        QPPv2y_te{ip,isd}{iy}=tez;
    end
end
end; clear T io ix tz z tez t

%% uncomment t-/ks-test whichever is to test & plot
pttsty_te=nan(nY,nY,3,nP,nSD,'single');
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:3
for iy1=1:nY
    te1=QPPv2y_te{ip,isd}{iy1}(:,ie);
    te1(isnan(te1))=[]; 
    for iy2=1:nY
        te2=QPPv2y_te{ip,isd}{iy2}(:,ie);
        te2(isnan(te2))=[]; 
        [h,p]=ttest2(te1,te2,'Alpha',alph,'vartype','unequal');
%         [h,p]=kstest2(te1,te2,'Alpha',alph);
        if ~h, pttsty_te(iy1,iy2,ie,ip,isd)=1;
        else, pttsty_te(iy1,iy2,ie,ip,isd)=p;
        end
    end
end
end
end
end; clear te1 te2 h p

s=[ss2(1) 1 ss2(3) ss(4)]; figure; set(gcf,'Position',s);
isd=1;
for ip=1:nP
    cnt=1;
for ie=[1 3] 
    p=pttsty_te(:,:,ie,ip,isd);
    subplot(nP,2,(ip-1)*2+cnt), im=imagesc(p); axis square; 
    set(im,'AlphaData',~isnan(p)); title(['QPP' num2str(ip)  ' t' nmte(cnt)]);
    xticks(1:nY); yticks(1:nY); xticklabels(YLB); yticklabels(YLB); 
    cnt=cnt+1;
end
end
saveas(gcf,[d2SAplt 'SAs_4_indep-ttst-RSN_SD' num2str(isd) '.png']); close
% saveas(gcf,[d2SAplt 'SAs_4_ks-tst-RSN_SD' num2str(isd) '.png']); close

