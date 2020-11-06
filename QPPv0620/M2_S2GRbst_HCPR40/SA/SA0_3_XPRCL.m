
%% Preprocessing: Comparing QPP1 obtained by analyzing parcels timeseries 
%% vs grayordinate timeseries  
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2bpr','p2qppf','p2u','p2prep','p2O','d2SA','d2SAplt',...
    'INFO','nsbj','nscn','nt','nT','nvx','nvxc','ivx','nVX',...
    'PL','PLh','nitr','ITP','ITPstp','tres','nX','PLe','PLc','tsh');
cth=[0.1 0.2]; ncth1=3; ITP=ITP{1}; stp=ITPstp; 
a={'QPP1Vx/','QPP1VxCrtx/'};
for i=1:2, a{i}=[d2SA a{i}]; if ~exist(a{i},'dir'), mkdir(a{i}); end; end
p2SV=cell(nsbj,2); for is=1:nsbj, for i=1:2
        p2SV{is,i}=[a{i} INFO.IDF{is,1}(1:6) '.mat']; end; end
p2s=[d2SA 'SA03.mat'];
load([p1 p2prep],'ixG');
QPP1p=load([p1 p2O],'QPPs'); QPP1p=QPP1p.QPPs(:,1); 
addpath([p1 p2qppf]); addpath(genpath([p1 p2u]));
set(0,'DefaultAxesTitleFontWeight','normal'); ss=get(0,'Screensize');

%% Detecting QPP1 in grayordinate
%%%% >> needs RAM & takes A LONG TIME! divide sbj across cpus/servers, also
%%%% consider using ITPl (limited # initial timepoints per scan); for
%%%% HCPR40, this analysis was done and results are available
% for is=1:nsbj    
%     D=zeros(nvx,nT,'single');
%     for iscn=1:nscn
%         load([p1 p2bpr{is,iscn}],'bpr');
%         D(:,(iscn-1)*nt+(1:nt))=zscore(bpr,[],2);
%     end; clear bpr
%     s=sprintf('Sbj%d',is);
%     QPPvx(D,nscn,PL,cth,ncth1,nitr,1,ITP,PLh,tres,p2SV{is,1},s,stp);
%     QPPvx(D(1:nvxc),nscn,PL,cth,ncth1,nitr,1,ITP,PLh,tres,p2SV{is,2},s,stp);
% end

%% Converting QPP1 from grayordinate to parcel space
QPP1v2p=cell(nsbj,2); QPP1v2p(:)={zeros(nX,PLe,'single')}; 
for i=1:2
for is=1:nsbj
    load(p2SV{is,i},'QPP');
    for ix=1:nX, QPP1v2p{is,i}(ix,:)=mean(QPP(ixG{ix},:)); end
end
end
save(p2s,'QPP1v2p');

%% Comparing QPP1s
c=zeros(nsbj,3);
for is=1:nsbj
    T=QPP1v2p(is,:);
    for i=1:2, c(is,i)=Tcomp(QPP1p{is},T{i},PLc,tsh); end
    c(is,3)=Tcomp(T{1},T{2},PLc,tsh);
end

cn=zeros(nsbj*(nsbj-1)/2,3); % comparing pairs of sbjs for null dist.
for i=1:3
    if i<=2, T=QPP1v2p(:,i); else, T=QPP1p; end
    ct=nan(nsbj);
for is1=1:nsbj-1    
for is2=is1+1:nsbj
    [ct(is1,is2),~]=Tcomp(T{is1},T{is2},PLc,tsh);         
end
end
    cn(:,i)=ct(~isnan(ct));
end
save(p2s,'c','cn','-append');

%% Plots
alm=[-0.5 0.5]; nsp1=5; nsp2=8; l={'Vx','VxCrtx','Prcl'};

for i=1:3
    if i<=2, T=QPP1v2p(:,i); else, T=QPP1p; end
    figure; set(gcf,'Position',ss);
    for is=1:nsbj, subplot(nsp1,nsp2,is), PLTT(T{is},PLc,alm); end
    saveas(gcf,[d2SAplt 'SA0_31_PRCL_effect_Analysis-' l{i} '.png']); close;
end

mc=round(myfshr(abs(c),1)*100)/100;
mcn=round(myfshr(abs(cn),1)*100)/100;
a={[l{1} ' vs ' l{3}],[l{2} ' vs ' l{3}],[l{1} ' vs ' l{2}]};
cbn=-1:0.01:1; figure; 
for i=1:3
    subplot(2,3,i), hist(c(:,i),cbn); title([a{i} ' ' num2str(mc(i))]);
    subplot(2,3,i+3), hist(cn(:,i),cbn); title(mcn(i)); 
end; subplot(2,3,1), title([a{1} ', med abs:' num2str(mc(1))]);
subplot(2,3,4), title(['Null dist, ' num2str(mcn(1))])
saveas(gcf,[d2SAplt 'SA0_32_PRCL_effect_absc_vs_abscn.png']); close

% figure; cnt=1; 
% for i=1:2
% for j=i+1:3
%     subplot(3,1,cnt), plot(c(:,i),c(:,j),'o'); 
%     axis square; axis([-1 1 -1 1]); grid on; cnt=cnt+1;
% end
% end
