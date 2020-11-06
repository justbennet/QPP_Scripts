
%% Detecting Group QPP by concatenating all scans of all subjects & ...
%% running the QPP algorithm for a limited random initial segments
%%
% The following steps are performed here:
% (1) Detecting the dominant QPP of the concatenated scans (QPPf1detect)
% (2) Phase-adjusting this QPP (QPPf2phadj) based on 1 or 2 seed parcels
% (3) Building QPP2 & above by averaging the contributing segments over the
% original scans NOT the residual scans (QPPf3xresid)
% (4) After detection of QPPi, regressing QPPs 1 to i scanwise(QPPf4regscn)
% to detect the next dominant QPP or QPPi+1 using the same steps
% > In the QPPf4regscn, FC after regression of QPPs 1-i is also calculated
% & later combined with FC of original scans (QPPf5dFC)
% (5) Building null patterns by averaging random segments (QPPf6null)

%% 
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nsbj','nscn','nscng','nX','nt','nTg','tres',...
    'p2B','p2O','p2Dr','p2qppf','nP','PL','cth','ncth1','nitr',...
    'ssg','ITPg','ITPgClk','ITPstp','PLh','tsh','PLc','cthph','sdph',...
    'nSD','ibY','iG2Y','fcbn','fcth','nnll');
addpath(p2qppf); ITPgc=[ITPg; ITPgClk];

%%
fprintf('Loading data\n'); 
load(p2B,'B');
D1=zeros(nX,nTg,'single');
for is=1:nsbj
	for i=1:nscn, D1(:,(is-1)*nt*nscn+(i-1)*nt+(1:nt))=B{is,i}; end 
end; clear B; D=D1;

QPP=cell(nP,1); TMX=cell(nP,1); C=zeros(nP,nTg,'single'); 
MET=zeros(nP,3,'single'); QPPa=cell(nP,nSD); TMXa=cell(nP,nSD);
METa=zeros(nP,nSD,3,'single'); cT1Tj=cell(nP,1); 
Cr=zeros(nP,nTg,'single'); FCr=cell(nP,1);

for ip=1:nP    
    
    tic; if ip~=1, D=Dr; clear Dr; end
    
    fprintf('QPP%d-f1detect\n',ip);
    [qpp,tmx,c,met,~,TMPL,TMXTMPL,CTMPL,SCMX]=QPPf1detect...
        (D,nscng,PL,cth{ip},ncth1(ip),nitr,ssg(ip),ITPg{ip},PLh,tres,...
        sprintf('QPP%d-f1detect',ip),ITPstp);
  
    fprintf('QPP%d-f2phadj\n',ip);
    [qppa,tmxa,~,meta]=QPPf2phadj...
        (qpp,TMPL,SCMX,TMXTMPL,CTMPL,tsh,PLc,cthph,sdph{ip},0,tres); 
    clear TMPL SCMX TMXTMPL CTMPL

    fprintf('QPP%d-f3xresid\n',ip);
    if ip~=1
        [qpp,c,met(1),qppa,~,meta(:,1)]=QPPf3xresid...
            (D1,tmx,tmxa,PLh,PLc,nscng,ssg(ip));
    end
    
    QPP{ip}=qpp; TMX{ip}=tmx; C(ip,:)=c; MET(ip,:)=met;
    QPPa(ip,:)=qppa; TMXa(ip,:)=tmxa; METa(ip,1:size(meta,1),:)=meta;
    
    fprintf('QPP%d-f4regscn\n',ip);
    [Dr,Cr(1:ip,:),FCr{ip}]=QPPf4regscn...
        (D1,QPP(1:ip),C(1:ip,:),nscng,PL,PLc,ibY,iG2Y,1);
    save(p2Dr{ip},'Dr','-v7.3'); 
    fprintf('QPP%d %dsec\n\n',ip,round(toc));
    
end; clear D Dr qpp tmx c met qppa tmxa meta
save(p2O,'QPP','TMX','C','MET','QPPa','TMXa','METa','Cr','ITPgc');

fprintf('QPPf5dFC\n');
[dFC,hFC,pFC,qFCN]=QPPf5dFC(D1,nscng,FCr,ibY,iG2Y,fcbn,fcth);
save(p2O,'dFC','hFC','pFC','qFCN','-append');

fprintf('QPPf6null\n'); 
[TMXN,qN]=QPPf6null(D1,MET(:,end),nscng,nnll,PL);
save(p2O,'TMXN','qN','-append');

