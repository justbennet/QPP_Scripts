
%% Plots in parcel-space
%%
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nscng','nt','p2O','p2qppf','p2u','d2plt','nP','PLc',...
    'PL','nY','G2Y','fcbn','fcth','nsd','nSD','nITPps','cthph'); 
load(p2O,'QPP','qN','dFC','hFC','pFC','qFCN','MET','QPPa','C','Cr');
addpath(p2qppf); addpath(p2u); % p2u needed for plot functions in p2qppf

alm=repmat([-1 1],nP,1); % default QPP amp. range when plotting
for i=4:nP, alm(i,:)=0.75*alm(i,:); end % set so only here
set(0,'DefaultAxesTitleFontWeight','normal');
ss=get(0,'Screensize'); ssw=ss; ssw(2)=ss(4)/4; ssw(4)=ss(4)/2;

%% Primary plots
%% 1 QPPs & FC within QPPs
figure; set(gcf,'Position',ssw);
for i=1:nP
    subplot(1,nP,i), PLTT(QPP{i},PLc,alm(i,:),1,1,~(i-1));
    title(['QPP' num2str(i)]); 
end; saveas(gcf,[d2plt '11_QPPs.png']); close;

figure; set(gcf,'Position',ssw);
for i=1:nP
    T=QPP{i}; T(sqrt(sum(T.^2,2))<qN(i,1),:)=nan; T(abs(T)<qN(i,2))=nan;
    subplot(1,nP,i), PLTT(T,PLc,alm(i,:),1,1,~(i-1));
    title(['QPP' num2str(i)]); 
end; saveas(gcf,[d2plt '11_QPPs_thrshld.png']); close;

figure; set(gcf,'Position',ssw);
for i=1:nP
    T=zeros(nY,PL,'single'); 
    for iy=1:nY, T(iy,:)=mean(QPP{i}(G2Y==iy,PLc),1); end; c=corr(T');
    c(abs(c)<=0.2)=nan;
    subplot(1,nP,i), PLTFCY7(c,[-1 1],0,1,~(i-1)); 
    title(['QPP' num2str(i)]);
end
a=get(gca,'Position'); a=[a(1)+a(3)+0.01 a(2)+a(4)/4 0.01 a(4)*0.5];
colorbar('position',a); saveas(gcf,[d2plt '12_FCwQPPs.png']); close

%% 2 FC reduction by regressing QPPs
fclm=[-0.6 0.8];
for i=1:nP
    figure; set(gcf,'Position',ss); PLTFC(dFC{i},fclm,1,1,1);
    saveas(gcf,[d2plt '21_dFC_QPPs1-' num2str(i) 'R.png']); close
end
c=dFC{nP}; c(abs(c)<qFCN)=nan;
figure; set(gcf,'Position',ss); PLTFC(c,fclm,1,1,1);
saveas(gcf,[d2plt '21_dFC_QPPs1-' num2str(nP) 'R_thrshld.png']); close

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

%% 4 QPPs basic metrics
a={'median C@maxima','median \Deltatmaxima (s)','# maxima'};
figure; set(gcf,'Position',ssw);
for i=1:3, subplot(1,3,i), plot(MET(:,i),'ko-'); 
    title(a{i}); xlabel('QPPs'); xlim([0 nP+1]); xticks(1:nP); end
saveas(gcf,[d2plt '41_Basic_Metrics.png']); close

%% Secondary plots
%% 5 QPP vs QPP-phase-adjusted
figure; set(gcf,'Position',ss);
for i=1:nP
    subplot(nSD+1,nP,i), PLTT(QPP{i},PLc,alm(i,:),1); 
    title(i); if ~(i-1), ylabel('QPP'); end
for isd=1:nsd(i)
if any(QPPa{i,isd})
    subplot(nSD+1,nP,isd*nP+i), PLTT(QPPa{i,isd},PLc,alm(i,:),1); 
    if ~(i-1), ylabel(['QPP phadj SD' num2str(isd)]); end
end
end
end; saveas(gcf,[d2plt '51_QPP_vs_QPPphadj.png']); close

%% 6 goodness of regression
cth1=0.2; cbn=cth1-0.05:0.01:0.8; h=zeros(nP,length(cbn),2);
for i=1:nP
    c=C(i,:); c=c(c>cth1); h(i,:,1)=hist(c,cbn);
    c=Cr(i,:); c=c(c>cth1); h(i,:,2)=hist(c,cbn);
end

a1={['Counts of correlation values above ' num2str(cth1)]; 
    'between QPP_i and original scans'}; 
a2={'~';'between QPP_i and residual scans'}; 
a=cell(1,nP); for i=1:nP, a{i}=['QPP' num2str(i)]; end
figure; set(gcf,'Position',ssw);
for i=1:2
    subplot(1,2,i), plot(cbn,h(:,:,i),'.-');
    title(eval(['a' num2str(i)])); xlabel('correlation'); legend(a);
    axis([cbn(1) cbn(end) 0 max(h(:))/i]);
end; saveas(gcf,[d2plt '61_Goodness_of_Regression.png']); close

figure; set(gcf,'Position',ss); nscn1=4; t=1:nscn1*nt;
for i=1:nP 
    subplot(nP,1,i), plot(C(i,t),'b'); hold on; plot(Cr(i,t),'r'); 
    title(['~ QPP' num2str(i)]); xlim([t(1) t(end)]);
    PLTUMET(nscn1*nt,nscn1,[-0.4 0.6],cth1);
end
subplot(nP,1,1), title('Correlation of QPP1');
legend('with original scan','with residual scan');
saveas(gcf,[d2plt '62_Goodness_of_Regression_2.png']); close
