
%% Obtaining group QPPs by combining subject QPPs
%%
% Across subjects, the following steps are performed here: 
% (1) All QPPs are combined (QPPscmbn)into a cell array & saved along a
% few cell arrays that contain phase-adjusted QPPs (QPPsa), timepoints of
% local maxima (TMXs & TMXsa) & basic metrics (METs & METsa); FCs after
% regression of QPPs 1-nP are averaged across sbj as a way of combining
% (2) FC before and after regression of QPPs are combined (QPPsdFC) 
% (3) Counts of transition between QPPs are calculated (QPPstrns), taking
% phase-adjusted QPPs with 1st seed and only 1st few QPPs (e.g., QPPs 1-3)
% (4) Phase-adjusted QPPs are averaged (QPPsfave) to obtain group QPP

%% 
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'p2S','p2B','p2O','p2G','p2qppf','ibY','iG2Y','fcbn','fcth',...
    'nT','nCNT','nP','nSD','nsd','nsbj','tsh','PLc','nitrg','cth','PLh',...
    'G2Y'); addpath(p2qppf);
fprintf('Loading Data'); tic; load(p2B,'B'); fprintf(' %ds\n',round(toc));

%%
fprintf('Combining'); tic;
[QPPs,TMXs,METs,QPPsa,TMXsa,METsa,FCr]=QPPscmbn(p2S);
fprintf(' %ds\n',round(toc)); 

fprintf('FC'); tic;
[dFC,hFC,pFC,qFCN,hFCN]=QPPsdFC(B,FCr,ibY,iG2Y,fcbn,fcth); 
fprintf(' %ds\n',round(toc)); 

fprintf('Transition Count'); tic;
[CNTtrns,CNTovlp,CNTnp]=QPPstrns(TMXsa(:,:,1),nT,nCNT);
fprintf(' %ds\n',round(toc)); 

c=cell(nP,nSD); QPP=c; TMX=c; CWQPP=c; clear c; 
NMX=zeros(nP,nSD,nsbj,'single');
for ip=1:nP
for isd=1:nsd(ip)
    fprintf('GrpQPP%dSD%d',ip,isd); tic;
    [QPP{ip,isd},TMX{ip,isd},NMX(ip,isd,:),CWQPP{ip,isd},...
        C,QPPg2s,Tg0,isref,cij,tshij,sgnij]=...
        QPPsfave(QPPsa(:,ip,isd),B,tsh,PLc,nitrg,cth{ip}(2),PLh,G2Y);
    save(p2G{ip,isd},'C','QPPg2s','Tg0','isref','cij','tshij','sgnij');
    fprintf(' %ds\n',round(toc));
end
end

save(p2O,'QPPs','TMXs','METs','QPPsa','TMXsa','METsa','dFC','hFC','pFC',...
    'qFCN','hFCN','CNTtrns','CNTovlp','CNTnp','QPP','TMX','NMX','CWQPP');
