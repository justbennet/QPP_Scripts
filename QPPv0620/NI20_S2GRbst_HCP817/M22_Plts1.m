
%% Primary plots in parcel-space
%%
clear; clc; close all; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nscn','nt','p2O','p2qppf','p2u','d2plt',...
    'nP','nSD','nsd','PLc','PL','nY','fcbn','fcth'); 
load(p2O,'QPP','CWQPP','dFC','hFC','pFC','stdFCN','METs','METsa',...
    'CNTtrns','CNTovlp','CNTnp');
addpath(p2qppf); addpath(p2u); % p2u needed for plot functions in p2qppf

alm=repmat([-1 1],nP,1); % QPP amp. range when plotting
for i=4:nP, alm(i,:)=0.5*alm(i,:); end
set(0,'DefaultAxesTitleFontWeight','normal','DefaultFigureColormap',jet);
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close
ssw=ss; ssw(2)=ss(4)/4; ssw(4)=ss(4)/2; 
sst=ss; sst(3)=round(sst(3)/2*1.1); fs=22;

%% 1 QPPs & FC within QPPs
figure; set(gcf,'Position',ss);
for i=1:nP
for isd=1:nsd(i)
    subplot(nSD,nP,(isd-1)*nP+i), PLTT(QPP{i,isd},PLc,alm(i,:),1,1,~(i-1));
    title(['QPP' num2str(i) 'SD' num2str(isd)]); 
end
end; saveas(gcf,[d2plt '11_QPPs.png']); close;

figure; set(gcf,'Position',ss);
for i=1:nP
for isd=1:nsd(i)
    subplot(nSD,nP,(isd-1)*nP+i), PLTFCY7(CWQPP{i,isd},[-1 1],0,1,~(i-1)); 
    title(['QPP' num2str(i) 'SD' num2str(isd)]);
end
end
subplot(nSD,nP,nP); a=get(gca,'Position'); 
a=[a(1)+a(3)+0.01 a(2)+a(4)/4 0.01 a(4)*0.5]; colorbar('position',a); 
saveas(gcf,[d2plt '12_FCwQPPs.png']); close

%% 2 FC reduction by regressing QPPs
fclm=[-0.6 0.8];
for i=1:nP
    figure; set(gcf,'Position',sst); PLTFC(dFC{i},fclm,1,1,1);
    set(gca,'fontsize',fs);
    f2s=[d2plt '21_dFC_QPPs1-' num2str(i) 'R.png'];
    print(gcf,f2s,'-dpng','-r1000'); close
end

hFC1=hFC./sum(hFC,2); hmx=max(hFC1(:)); 
xtk=[fclm(1) -fcth 0 fcth fclm(2)]; ix4p={1:2,4:5};
cl=[0 0.45 0.65 0.85]; clp='br'; 
l={'QPPs 1-3 regressed','QPPs 1-2 regressed',...
    'QPP1 regressed'};
s=[ss2(1) 1 1.5*ss2(3) ss(4)*0.7];
figure; set(gcf,'Position',s); hold on
for i=nP+1:-1:1
    plot(fcbn,hFC1(i,:),'color',(i-1)/(nP+1)*ones(1,3),'linewidth',4);
end
for j=xtk(2:end-1), plot([j j],[0 hmx],'k'); end
axis([fclm 0 hmx]); xticks(xtk); yticks([]); box on; 
xlabel('correlation'); set(gca,'fontsize',fs+4);
for i=1:2
    a=xtk(ix4p{i});
    i1=find(abs(fcbn-a(1))<0.001); i2=find(abs(fcbn-a(2))<0.001);
    patch([a(1):0.01:a(2) a(2) a(1)],[hFC1(1,i1:i2) 0 0],clp(i),'facealpha',0.3);
end
legend(l);
print(gcf,[d2plt '22_histFC_4ppr.png'],'-dpng','-r300'); close

clr='brk';
l={'\Delta(c<-','\Delta(c>+','\Delta(|c|>'};
for i=1:3, l{i}=[l{i} num2str(fcth) ')']; end
s=[ss2(1) 1 ss2(3)*0.8 ss(4)]; 
figure; set(gcf,'Position',s); hold on; 
for i=1:3, plot(pFC(:,i),[clr(i) 'o-'],'linewidth',4,'markersize',10); end
xlim([0 nP+2]); xticks(1:nP+1); xticklabels(0:nP); xlabel('regressed QPPs'); 
grid on; box on; set(gca,'fontsize',fs+8); legend(l); 
print(gcf,[d2plt '23_pdFC_4ppr.png'],'-dpng','-r300'); close

pdFC=round((pFC(1,:)-pFC(end,:))./pFC(1,:)*100);

%% 3 QPPs in grayordinate & time-of-peak maps: plotted by other Plts script

%% 4.1 QPPs basic Metrics
mMET=squeeze(mean(METs)); sMET=squeeze(std(METs));
mMETa=squeeze(mean(METsa)); sMETa=squeeze(std(METsa));

a={'median C@maxima','median \Deltatmaxima (s)','# maxima'};
figure; set(gcf,'Position',ssw);
for i=1:3, subplot(1,3,i), hold on;
    for isd=1:nsd(i), errorbar(mMETa(i,:,isd),sMETa(i,:,isd),'o-'); end
    errorbar(mMET(i,:),sMET(i,:),'ko-');
    title(a{i}); xlabel('QPPs'); xlim([0 nP+1]); xticks(1:nP); box on
end
saveas(gcf,[d2plt '41_Basic_Metrics.png']); close

%% 4.2 Transition count between QPPs
% % CNTtrns: transition count from QPPi in rows to QPPj in columns
% % CNTnp: #QPPs per included sbjs (i.e., sbjs having all phase-adj. QPPs)
% % CNTovlp: overlap count between QPPs i&j
% d1=sum(CNTnp)-sum(CNTtrns);
% % note overlapping QPPs & last QPP per sbj aren't counted (see next)
% d2=sum(d1)-size(CNTnp,1)-sum(CNTovlp(:)); % zero as expected

Ptrns=round(CNTtrns./sum(CNTtrns,2)*100)/100;
X={CNTtrns,Ptrns};
figure; 
for i=1:2
    subplot(1,2,i), imagesc(X{i}); colormap('jet'); colorbar; 
    xticks(1:nP); yticks(1:nP); axis square; 
end
saveas(gcf,[d2plt '42_TrnsCnt_TrnsProb.png']); close
