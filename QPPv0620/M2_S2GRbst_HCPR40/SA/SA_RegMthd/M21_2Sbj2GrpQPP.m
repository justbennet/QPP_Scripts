
%%
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'p2S','p2B','p2O','p2qppf','ibY','iG2Y','fcbn','fcth',...
    'nT','nP','nSD','nsd','nsbj','tsh','PLc','nitrg','cth','PLh',...
    'G2Y','nnll','PL'); addpath(p2qppf);
fprintf('Loading Data'); tic; load(p2B,'B'); fprintf(' %ds\n',round(toc));

%% 
fprintf('Combining'); tic;
[QPPs,TMXs,METs,QPPsa,TMXsa,METsa,FCr]=QPPscmbn(p2S);
fprintf(' %ds\n',round(toc)); 

fprintf('FC'); tic;
[dFC,hFC,pFC,qFCN,hFCN]=QPPsdFC(B,FCr,ibY,iG2Y,fcbn,fcth); 
fprintf(' %ds\n',round(toc)); 

c=cell(nP,nSD); QPP=c; TMX=c; CWQPP=c; clear c; 
NMX=zeros(nP,nSD,nsbj,'single');
for ip=1:nP
for isd=1:nsd(ip)
    fprintf('GrpQPP%dSD%d',ip,isd); tic;
    [QPP{ip,isd},TMX{ip,isd},NMX(ip,isd,:),CWQPP{ip,isd}]=...
        QPPsfave(QPPsa(:,ip,isd),B,tsh,PLc,nitrg,cth{ip}(2),PLh,G2Y);
    fprintf(' %ds\n',round(toc));
end
end

fprintf('Null Patterns\n');
[TMXN,qN]=QPPsnull(B,max(sum(NMX,3),[],2),nnll,PL);

save(p2O,'QPPs','TMXs','METs','QPPsa','TMXsa','METsa','dFC','hFC','pFC',...
    'qFCN','QPP','TMX','NMX','CWQPP','TMXN','qN');
