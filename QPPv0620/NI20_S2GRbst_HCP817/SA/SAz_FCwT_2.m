
%% Supp analyses for FC within QPPs - part2
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2O','nX','nP','p2qppf','PLc','p2u','nY','ibY','iG2Y',...
    'd2SA','p2B','nsbj','nscn','nt','nT','PLe','PL','PLh','p2S','d2SAplt'); 
p2O=[p1 p2O]; p2B=[p1 p2B]; p2qppf=[p1 p2qppf]; p2u=[p1 p2u];
addpath(p2qppf); addpath(genpath(p2u));
ss=get(0,'Screensize');

iu=find(triu(ones(nX),1)); nX2=nX*(nX-1)/2;

%% Null value for FC within QPPs
for k=1:2
    
if k==1, load([d2SA 'SA2.mat'],'PNP'); % averages of random segments
elseif k==2, load([d2SA 'SAz_FCwT_2.mat'],'PNP'); % see next section
end

[nrep,n]=size(PNP);
FCWTN=nan(nX,nX,nrep,n); fcwtn=zeros(nX2,nrep,n);
for ir=1:nrep
for in=1:n
    T=PNP{ir,in}(:,PLc);
    To=T; for i=1:nY, To(ibY(i)+1:ibY(i+1),:)=T(iG2Y{i},:); end
    c=corr(To'); FCWTN(:,:,ir,in)=c;
    fcwtn(:,ir,in)=c(iu);
end
end

nrep1=16; a=randperm(nrep); a=sort(a(1:nrep1)); 
figure; set(gcf,'Position',ss); cnt=1; 
for in=1:n
for ir=a
    subplot(n*2,nrep1/2,cnt), PLTFC(FCWTN(:,:,ir,in),[-1 1],0,0,0); cnt=cnt+1;
end
end

cn=reshape(fcwtn,nX2,nrep*n);
cn=corr(cn); iu1=find(triu(ones(nrep*n),1)); cn=cn(iu1);
cnst=[min(cn) median(cn) std(cn) max(cn)];
figure; hist(cn,0:0.01:0.5);

end

%% Phase-randomized timeseries, averaging random segments & QPPs' segments
% load(p2B,'B');
% BN=cell(nsbj,nscn);
% for is=1:nsbj
%     is
% for iscn=1:nscn
%     amp=abs(fft(B{is,iscn}')); 
%     r=rand(nt,nX,'single'); r=zscore(r); pha=angle(fft(r));  
%     N=ifft(amp.*exp(sqrt(-1)*pha)); BN{is,iscn}=real(N)';
% end
% end
% B=BN; clear BN; 
% save([d2SA 'SAz_BN.mat'],'B','-v7.3');

% load([d2SA 'SA2.mat'],'TN','nmxn'); 
% addpath('./SA2_SigActiv/');
% PNP=fAmpNullPrcl(TN,nmxn,[d2SA 'SAz_BN.mat']);
% save([d2SA 'SAz_FCwT_2.mat'],'PNP');
% 
% load(p2O,'TMX','NMX'); 
% QPPN=cell(nP,1); QPPN(:)={zeros(nX,PLe,'single')};
% for is=1:nsbj
%     is
%     D=zeros(nX,nT,'single'); 
%     for i=1:nscn, D(:,(i-1)*nt+(1:nt))=B{is,i}; end
% for ip=1:nP
%     tmx=TMX{ip,1}{is};
%     if any(tmx), T=Tbld(D,tmx,PL,PLh,0); QPPN{ip}=QPPN{ip}+T; end
% end
% end
% for ip=1:nP, QPPN{ip}=QPPN{ip}/sum(NMX(ip,1,:)); end
% save([d2SA 'SAz_FCwT_2.mat'],'QPPN','-append');

