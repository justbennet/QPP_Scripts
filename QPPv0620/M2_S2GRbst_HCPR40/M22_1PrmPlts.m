
%% Primary plots in parcel-space
%%
clear; clc; close all; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nscn','nt','p2O','p2qppf','p2u','d2plt',...
    'nP','nSD','nsd','PLc','PL','nY','fcbn','fcth'); 
load(p2O,'QPP','CWQPP','dFC','hFC','pFC','qFCN','METs','METsa',...
    'CNTtrns','CNTovlp','CNTnp');
addpath(p2qppf); addpath(p2u); % p2u needed for plot functions in p2qppf

alm=repmat([-1 1],nP,1); % QPP amp. range when plotting
for i=4:nP, alm(i,:)=0.5*alm(i,:); end
set(0,'DefaultAxesTitleFontWeight','normal','DefaultFigureColormap',jet);
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close
ssw=ss; ssw(2)=ss(4)/4; ssw(4)=ss(4)/2; 
sst=ss; sst(3)=round(sst(3)/2*1.1);

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
for i=[1 nP]
    figure; set(gcf,'Position',ss); PLTFC(dFC{i},fclm,1,1,1);
    saveas(gcf,[d2plt '21_dFC_QPPs1-' num2str(i) 'R.png']); close
end

figure; plot(fcbn,hFC./sum(hFC,2),'.-');
title('Histogram of FC values'); xlabel('correlation'); yticklabels([]);
xlim(fclm); xticks([fclm(1) -fcth 0 fcth fclm(2)]); grid on; box on; 
a={'','QPP1 regressed'}; 
for i=3:nP+1, a{i}=['QPPs 1-' num2str(i-1) ' regressed']; end; legend(a);
saveas(gcf,[d2plt '22_histFC.png']); close

figure; clr='brk'; hold on; 
for i=1:3, plot(round(pFC(:,i)),[clr(i) 'o-']); end
title(['Percentage of FC values with magnitude above ' num2str(fcth)]);
xlabel('regressed QPPs'); axis square;
xlim([0 nP+2]); xticks(1:nP+1); xticklabels(0:nP); grid on; box on; 

pdFC=round((pFC(1,:)-pFC(end,:))./pFC(1,:)*100);
a={'\Delta(FC<-','\Delta(FC>+','\Delta(|FC|>'};
for i=1:3, a{i}=[a{i} num2str(fcth) ')=' num2str(pdFC(i)) '%']; end
legend(a); saveas(gcf,[d2plt '23_pdFC.png']); close

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
