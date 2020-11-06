
%% Significant activation/deactivation within QPP
%% Part2N) Finding magnitude of null patterns (99th quantile can sever as
%% a magnitude threhold for QPPs' timecourses)
%%
clear; clc; close all; p2s='SA2.mat'; p2='../../'; p1='../';
p2p=dir([p2 'Params_*.mat']); p2p=[p2 p2p.name];
load(p2p,'d2SA','nvx','PL','PLc','d2SAplt','irgn','nrgn','nmrgn1');
d2SA=[p1 d2SA(3:end)]; d2SAplt=[p1 d2SAplt(3:end)];
load([d2SA p2s],'PN'); [nrep,n]=size(PN);
cl='gk'; lw=1.5; fs=20; ss=get(0,'Screensize'); 
set(0,'DefaultAxesTitleFontWeight','normal');

%%
AG=cell(nrgn,n); AN=zeros(nvx*PL*nrep,n,'single'); 
for in=1:n
    a=zeros(nvx,PL,nrep,'single'); 
    for ir=1:nrep, a(:,:,ir)=PN{ir,in}(:,PLc); end
    for ig=1:nrgn, z=a(irgn{ig},:,:); AG{ig,in}=z(:); end
    AN(:,in)=a(:);
end; clear a z
amx=ceil(max(abs(AN(:)))*100)/100;
abn=-amx:0.001:amx; lbn=length(abn);

%%
h=zeros(lbn,nrgn,n); q=zeros(nrgn,n); 
h1=zeros(lbn,n); q1=zeros(nrgn,1); 
for in=1:n
    for ig=1:nrgn
        z=AG{ig,in}; h(:,ig,in)=hist(z,abn); q(ig,in)=quantile(abs(z),0.99); 
    end
    z=AN(:,in); h1(:,in)=hist(z,abn); q1(in)=quantile(abs(z),0.99); 
end

%%
figure; s=get(gcf,'Position'); s1=s; s1(3)=1.5*s1(3); 
set(gcf,'Position',s1); hold on; hmx1=max(h1(:));
for i=1:n
    plot(abn,h1(:,i),cl(i),'linewidth',lw); 
    for j=[-1 1]
        z=j*q1(i); plot([z z],[0 hmx1],[cl(i) '--'],'Linewidth',lw);
    end
end; axis([-amx amx 0 hmx1]); box on; xticks(-0.05:0.025:0.05)
ylabel('count of vertices/voxels'); 
xlabel('amplitude of ~20 timecourse per vertex/voxel');
legend('N1 random segments','99% quantile','',...
    'N3 random segments','99% quantile',''); set(gca,'fontsize',fs); 
saveas(gcf,[d2SAplt 'SA2_2n_1_AmpThrsh.png']); close

figure; s2=s; s2([2 4])=ss([2 4]); set(gcf,'Position',s2);
for ig=1:nrgn
    hf=h(:,ig,:); hmx=max(hf(:));
    subplot(nrgn,1,ig),hold on
for i=1:n
    plot(abn,h(:,ig,i),cl(i),'linewidth',lw); 
    for j=[-1 1]
        z=j*q(ig,i); plot([z z],[0 hmx],[cl(i) '--'],'Linewidth',lw);
    end
    axis([-amx amx 0 hmx]); ylabel(nmrgn1{ig}); xticks(-0.05:0.025:0.05)
    if ig~=nrgn, xticklabels([]); end; yticklabels([]);
    set(gca,'fontsize',fs); box on;  
end
end; xlabel('amplitude of ~20s timecourse');
saveas(gcf,[d2SAplt 'SA2_2n_2_AmpTrsh.png']); close
