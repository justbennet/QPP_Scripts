
%% Preprocessing: averaging over GM&WM&CSF-regressed scans (QPPs) vs 
%% WM&CSF-regressed scans (i.e., effect of GSR)
%%
clear; clc; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name]; 
load(p2p,'p2prep','p2O','p2qppf','p2u','d2SAplt','nsbj',...
    'nscn','nX','nt','nT','nP','nSD','nsd','PL','PLh','PLe','PLc','tsh');
load([p1 p2prep],'p2BWCR'); a=strfind(p2prep,'/'); 
p2B=[p1 p2prep(1:a(end)) p2BWCR]; load(p2B,'B');
load([p1 p2O],'QPPs','TMXs','QPPsa','TMXsa'); 
addpath([p1 p2qppf]); addpath([p1 p2u]);
set(0,'DefaultAxesTitleFontWeight','normal');

%%
c=nan(nP,nSD+1,nsbj);
for is=1:nsbj
    D=zeros(nX,nT,'single');
    for i=1:nscn, D(:,(i-1)*nt+(1:nt))=B{is,i}; end
for ip=1:nP
    T=Tbld(D,TMXs{is,ip},PL,PLh,1);
% 	c(ip,1,is)=Tcomp(QPPs{is,ip},T,PLc,tsh); % using more strict corr instd
    P=QPPs{is,ip}(:,PLc); T=T(:,PLc); c(ip,1,is)=corr(P(:),T(:));
for isd=1:nsd(ip)
    tmx=TMXsa{is,ip,isd}; n=length(tmx);
    if n
        T=Tbld(D,tmx,PL,PLh,1); 
        % c(ip,isd+1,is)=Tcomp(QPPsa{is,ip,isd},T,PLc,tsh); 
        P=QPPsa{is,ip,isd}(:,PLc); T=T(:,PLc); c(ip,isd+1,is)=corr(P(:),T(:));
    end
end
end
end
[mc,~,mcn,mcx]=myfshr(c,3); [~,n]=myst(c(:),2);

figure; 
subplot(1,2,1), hist(c(:),n:0.01:1); 
title({'QPP built by ave GM&WM&CSF-reg. scans vs'; 
    '~ WM&CSF-reg. scans, all QPPs all sbj'});
xlabel('corr'); xlim([n-0.01 1.01]); axis square; 
subplot(1,2,2), PLTC(mc,[mcn mcx],1,{'','phadjSD1','~2'},1:nP); 
title('med corr btwn pairs'); ylabel('QPP');
saveas(gcf,[d2SAplt 'SA0_1_GMR-effect.png']); close

