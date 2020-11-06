
%% Supp analyses for FC change by regressing QPPs & FC within QPPs
%%
clear; clc; close all; p1='../';
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2O','nX','nP','p2qppf','PLc','p2u','nY','ibY','iG2Y',...
    'd2SA','p2B','nsbj','nscn','nt','nT','PLe','PL','PLh','p2S','d2SAplt'); 
p2O=[p1 p2O]; p2B=[p1 p2B]; p2qppf=[p1 p2qppf]; p2u=[p1 p2u];
addpath(p2qppf); addpath(genpath(p2u));
set(0,'DefaultAxesTitleFontWeight','normal'); fs=20;
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close;

il=find(tril(ones(nX),-1)); iu=find(triu(ones(nX),1)); nX2=nX*(nX-1)/2;

%% FC matrices: Variance & correlation btwn pairs
%%
load(p2O,'dFC'); 
fc=zeros(nX2,nP+1,'single');
fc(:,1)=dFC{1}(iu); for i=1:nP, fc(:,i+1)=dFC{i}(il); end

fcv=round(std(fc).^2,3);
fcvd=1-fcv/fcv(1);

[fcc,fccp]=corr(fc); 

%% Variance of FC matrices using shuffled correlation timecourse of QPPs
%%
% load(p2B,'B');
% FCr=cell(nP,1); FCr(:)={zeros(nX)};
% for is=1:nsbj
%     is
%     D=zeros(nX,nT,'single');
%     for iscn=1:nscn, D(:,(iscn-1)*nt+(1:nt))=B{is,iscn}; end
%     TT=cell(nP,1); CT=zeros(nP,nT,'single');
%     for ip=1:nP
%         load([p1 p2S{is,ip}],'QPP','C'); TT{ip}=QPP; 
%         CT(ip,:)=C(randperm(nT));
%         [~,~,fcr]=QPPf4regscn(D,TT(1:ip),CT(1:ip,:),nscn,PL,PLc,ibY,iG2Y,0);
%         FCr{ip}=FCr{ip}+fcr;
%     end
% end
% for ip=1:nP, c=FCr{ip}/nsbj; c=(exp(2*c)-1)./(exp(2*c)+1); c(c>=0.9999)=1; 
%     FCr{ip}=c; end
% save([d2SA 'SAs_dFC_ShffC.mat'],'FCr');

load([d2SA 'SAs_dFC_ShffC.mat'],'FCr');
fcsh=zeros(nX2,nP+1,'single');
fcsh(:,1)=dFC{1}(iu); for i=1:nP, fcsh(:,i+1)=FCr{i}(il); end
fcshv=std(fcsh).^2;
fcshvd=1-fcshv/fcshv(1);

%% FC within QPPs
%%
load(p2O,'QPP');
FCWT=nan(nX,nX,nP); fcwt=zeros(nX2,nP);
for ip=1:nP
    T=QPP{ip,1}(:,PLc);
    To=T; for i=1:nY, To(ibY(i)+1:ibY(i+1),:)=T(iG2Y{i},:); end
    c=corr(To'); FCWT(:,:,ip)=c;
    fcwt(:,ip)=c(iu);
end

figure; s=[1 ss2(2) ss(3) ss2(4)]; set(gcf,'Position',s);
for ip=1:nP 
    subplot(1,nP,ip), PLTFC(FCWT(:,:,ip),[-1 1],0,1,~(ip-1)); 
    set(gca,'fontsize',fs); xtickangle(90);
end
f2s=[d2SAplt 'SAs_6_FCwQPP.png'];
print(gcf,f2s,'-dpng','-r200'); close

FCWTc=corr(fcwt);
