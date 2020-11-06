
%% Significant activation/deactivation within QPP
%% Part3N) Magnitude threshold for parcels/networks
%%
clear; clc; close all; p2s='SA2.mat'; p2='../../'; p1='../';
p2p=dir([p2 'Params_*.mat']); p2p=[p2 p2p.name];
load(p2p,'d2SA','p2u','ivx','nY','PL','irf','PLc','d2SAplt'); 
d2SA=[p1 d2SA(3:end)]; d2SAplt=[p1 d2SAplt(3:end)];
load([d2SA p2s],'PN'); [nrep,n]=size(PN); p2u=[p2 p2u];
load([p2u 'myGlssr.mat'],'Y','nG','ixG'); Y=Y(ivx); 
nX=1+nG; iX=[irf{1}; ixG];
ss=get(0,'Screensize'); cl='gk'; lw=1.5; fs=20; 
set(0,'DefaultAxesTitleFontWeight','normal');

%%
y=zeros(nY,PL,nrep,n,'single'); 
rf=zeros(nX,PL,nrep,n,'single'); 
for in=1:n
for ir=1:nrep
    z=zeros(nY,PL,'single');
    for i=1:nY, z(i,:)=mean(PN{ir,in}(Y==i,PLc)); end 
    y(:,:,ir,in)=z;
    z=zeros(nX,PL,'single');
    for i=1:nX, z(i,:)=mean(PN{ir,in}(iX{i},PLc)); end
    rf(:,:,ir,in)=z;
end
end

ayn=reshape(y,[],n);
amx=max(abs(ayn(:))); abn=-amx:0.0001:amx; lbn=length(abn);

h=zeros(lbn,n); for i=1:n, h(:,i)=hist(ayn(:,i),abn); end; hmx=max(h(:)); 
qy=quantile(ayn,0.99);

figure; hold on; 
for i=1:n
    plot(abn,h(:,i),cl(i),'linewidth',lw); 
    for j=[-1 1]
        z=j*qy(i); plot([z z],[0 hmx],[cl(i) '--'],'Linewidth',lw);
    end
end; axis([-amx amx 0 hmx]); box on; set(gca,'fontsize',fs);
xlabel('amplitude of ~20s timecourse per RSN'); ylabel('Count of RSNs')
saveas(gcf,[d2SAplt 'SA2_3n_1_AmpThsh-RSN.png']); close

%%
rfn=reshape(rf,[],n);
rmx=max(abs(rfn(:))); rbn=-rmx:0.001:rmx; lbn=length(rbn);

h=zeros(lbn,n); for i=1:n, h(:,i)=hist(rfn(:,i),rbn); end; hmx=max(h(:));
qr=quantile(rfn,0.99);

figure; hold on; 
for i=1:n
    plot(rbn,h(:,i),cl(i),'linewidth',lw); 
    for j=[-1 1]
        z=j*qr(i); plot([z z],[0 hmx],[cl(i) '--'],'Linewidth',lw);
    end
end; axis([-rmx rmx 0 hmx]); box on; set(gca,'fontsize',fs);
xlabel('amplitude of ~20s timecourse per parcel'); ylabel('Count of parcels')
saveas(gcf,[d2SAplt 'SA2_3n_1_AmpThsh-Parcel.png']); close

%%
