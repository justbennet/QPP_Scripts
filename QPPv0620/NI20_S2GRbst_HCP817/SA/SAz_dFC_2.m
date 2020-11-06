
%% Supp analyses for FC change by regressing QPPs - part2
%% Change in FC matrices within blocks corresponding to pairs of RSNs
%%
clear; clc; close all; p1='../';
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'nX','nY','ibY','nP','p2O','YLB'); p2O=[p1 p2O];
load(p2O,'dFC'); 

%% Finding indices of blocks corresponding to RSN-pairs
nY2=nY*(nY-1)/2+nY;
z=zeros(nX,nX,nY2); I=zeros(nY2,2); cnt=1; 
for iy1=1:nY
for iy2=iy1:nY
    z(ibY(iy1)+1:ibY(iy1+1),ibY(iy2)+1:ibY(iy2+1),cnt)=1;
    I(cnt,:)=[iy1 iy2];
    cnt=cnt+1;
end
end; clear iy1 iy2
for i=1:nX, for j=1:nX, z(j,i,:)=z(i,j,:); end; end

il2=cell(nY2,1); iu2=cell(nY2,1);
for i=1:nY2
    il2{i}=find(tril(z(:,:,i),-1)); iu2{i}=find(triu(z(:,:,i),1)); 
end; clear z

%% Breaking FC matrices into blocks of RSN-pairs
fc2=cell(nY2,1); 
for i=1:nY2, fc2{i}(:,1)=dFC{1}(iu2{i}); end
for ip=1:nP, for iy=1:nY2, fc2{iy}(:,1+ip)=dFC{ip}(il2{iy}); end; end

%% Varience of FC matrices per block
fc2v=zeros(nY2,nP+1);
for i=1:nY2, fc2v(i,:)=std(fc2{i}).^2; 
    fc2v(i,:)=round(fc2v(i,:)/fc2v(i,1)*100,1); end

[~,is]=sort(fc2v(:,end)); % sorting 28 blocks based on drop in variance
J=I(is,:);
K=nan(nY); for i=1:nY2, K(J(i,2),J(i,1))=nY2-i+1; end
K(5,:)=nan; K(:,5)=nan;

figure; im=imagesc(K); % sorted order of blocks (1:28) low->high var-drop
axis square; colormap('jet'); colorbar
set(im,'AlphaData',~isnan(K)); xticks(1:nY); yticks(1:nY);
xticklabels(YLB); yticklabels(YLB);

%% Correlation of FC-blocks after regressing QPPs
c2=cell(nY2,1); id=find(eye(nP+1));
for i=1:nY2, c=corr(fc2{i}); c(id)=nan; c2{i}=c; end; clear c

ss=get(0,'screensize'); s=ss; s(3)=s(4);
figure; set(gcf,'Position',s)
for i=1:nY2
    c=c2{i};
    subplot(nY,nY,(I(i,1)-1)*nY+I(i,2)), im=imagesc(c,[-1 1]); axis square;
    colormap('jet'); set(im,'AlphaData',~isnan(c));
    xticklabels([]); yticklabels([]); 
end
