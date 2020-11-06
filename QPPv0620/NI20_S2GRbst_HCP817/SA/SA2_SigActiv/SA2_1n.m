
%% Significant activation/deactivation within QPP
%% Part1N) Building null patterns & finding root sum square (~power) of
%% timecourses (99th quantile can be power threhold for QPPs' timecourses)
%%
clear; clc; close all; nIS=4; p2s='SA2.mat'; p2='../../'; p1='../';
p2p=dir([p2 'Params_*.mat']); p2p=[p2 p2p.name]; 
load(p2p,'nsbj','d2SA','p2B','PLc','PL','d2SAplt','nrgn','irgn','nmrgn1'); 
d2SA=[p1 d2SA(3:end)]; d2SAplt=[p1 d2SAplt(3:end)]; p2B=[p2 p2B];

IS=cell(nIS,1); p2IS=IS; r=floor(nsbj/nIS);
for i=1:nIS, IS{i}=(1:r)+(i-1)*r; p2IS{i}=[d2SA 'f' num2str(i) '_' p2s]; end
IS{nIS}=[IS{nIS} IS{nIS}(end)+1:nsbj];

load([d2SA p2s],'TN','nmxn'); [nrep,n]=size(TN);
ss=get(0,'screensize'); cl='gk'; lw=1.5; fs=20;

%% Building null patterns by averaging random segments
%% fAmpNull: Run on servers (i.e., nIS nodes for i=1:nIS)
% i=1; fAmpNull(IS{i},p2IS{i},TN);

%% Combining 
% load(p2IS{1},'PN'); 
% for i=2:nIS
%     a=load(p2IS{i},'PN');
%     for ir=1:nrep, for in=1:n 
%             PN{ir,in}=PN{ir,in}+a.PN{ir,in}; end; end
% end
% for ir=1:nrep, for in=1:n, PN{ir,in}=PN{ir,in}/nmxn(in); end; end
% save([d2SA p2s],'PN','-v7.3');
% save([d2SA p2s],'TN','nmxn','-append');

%% Also averaging in parcel-space
% PNP=fAmpNullPrcl(TN,nmxn,p2B);
% save([d2SA p2s],'PNP','-append');

%% Finding root sum square of timecourses of the null patterns
%%
load([d2SA p2s],'PN');
load(p2p,'nvx','nX','nvxc');

pn=zeros(nvx,nrep,n,'single'); 
for ir=1:nrep
for in=1:n
    pn(:,ir,in)=sqrt(sum(PN{ir,in}(:,PLc).^2,2)); 
end
end

%% Finding stats to serve as power threhold for QPPs
%%
p=pn; 
pmx=ceil(max(p(:))*100)/100; pbn=0:0.001:pmx; lbn=length(pbn);
q=zeros(nrgn,n); h=zeros(lbn,nrgn,n);
q1=zeros(1,n); h1=zeros(lbn,n);
for i=1:n
    for ig=1:nrgn
        a=p(irgn{ig},:,i); a=a(:); 
        q(ig,i)=quantile(a,0.99); h(:,ig,i)=hist(a,pbn);
    end
    a=p(:,:,i); a=a(:); q1(i)=quantile(a,0.99); h1(:,i)=hist(a,pbn);
end
hmx1=max(h1(:)); 

%% 
figure; s=get(gcf,'Position'); s1=s; s1(3)=1.5*s1(3); 
set(gcf,'Position',s1); hold on; 
for i=1:n
    plot(pbn,h1(:,i),cl(i),'Linewidth',lw); 
    plot([q1(i) q1(i)],[0 hmx1],[cl(i) '--'],'Linewidth',lw);
end
ylabel('count of vertices/voxels'); 
xlabel('root sum square of ~20s-long timecourse');
legend('N1 random segments','99% quantile',...
    'N3 random segments','99% quantile');
axis([0 pbn(end) 0 hmx1]); set(gca,'fontsize',fs); box on
saveas(gcf,[d2SAplt 'SA2_1n_1_PwrThrsh.png']); close

figure; s2=s; s2([2 4])=ss([2 4]); set(gcf,'Position',s2);
for ig=1:nrgn
    hf=h(:,ig,:); hmx=max(hf(:));
    subplot(nrgn,1,ig), hold on 
    for i=1:n
        plot(pbn,h(:,ig,i),cl(i),'Linewidth',lw); 
        plot([q(ig,i) q(ig,i)],[0 hmx],[cl(i) '--'],'Linewidth',lw);
    end
    axis([0 pbn(end) 0 hmx]); yticklabels([]); ylabel(nmrgn1{ig}); 
    if ig~=nrgn, xticklabels([]); end; set(gca,'fontsize',fs); box on
end
xlabel('root sum square of ~20s timecourse');
saveas(gcf,[d2SAplt 'SA2_1n_2_PwrThrsh-rgn.png']); close
