clear; clc

%%
addpath('./cifti-matlab-master');
e=ft_read_cifti('empty.dtseries.nii');
es=e.brainstructure; e1=e.dtseries; en=e.brainstructurelabel';

nVX=size(es,1); 
ivx=single(find(~isnan(e1))); 
nvx=length(ivx);
nVXC=find(es==3); nVXC=nVXC(1)-1;
es1=es(ivx); nvxc=find(es1==3); nvxc=nvxc(1)-1;
ivxc=ivx(1:nvxc);
save('myHCPcft.mat','nVX','ivx','nvx','nVXC','ivxc','nvxc');

%%
sc{1}=[10 11];   nmsc{1}='Cerebellum';     nmsc1{1}='C'; 
sc{2}=[20 21];   nmsc{2}='Thalamus';       nmsc1{2}='T';
sc{3}=[14 15];   nmsc{3}='Hippocampus';    nmsc1{3}='H';
sc{4}=[5 6];     nmsc{4}='Amygdala';       nmsc1{4}='A';
sc{5}=[7 12 13 16 17]; nmsc{5}='Brainstem & Deep brain';  nmsc1{5}='BD';  
sc{6}=[8 9 18 19 3 4]; nmsc{6}='Striatum'; nmsc1{6}='S';
sc=sc'; nmsc=nmsc'; nmsc1=nmsc1';

ep1=e.pos(ivx,:); 
N=length(sc); ivxsc=cell(N,1); pvxsc=cell(N,1); 
for i=1:N
    ind=[]; 
    for j=1:length(sc{i}), ind=[ind; single(find(es1==sc{i}(j)))]; end
    ivxsc{i}=ind; pvxsc{i}=single(ep1(ind,:));
end
save('myHCPcft.mat','sc','nmsc','nmsc1','ivxsc','pvxsc','-append')

% for i=1:N % for checking
%     m=nan(nvx,1); m(ivxsc{i})=1; M=nan(nVX,1); M(ivx)=m; e.dtseries=M;
%     ft_write_cifti(['./U3myPrcls/HCPMSK_' num2str(i) nmsc1{i}],e,...
%         'parameter','dtseries');
% end

nmrgn=[{'Cortex'}; nmsc]; nmrgn1=[{'Cx'}; nmsc1]; 
irgn=[single(1:nvxc)'; ivxsc]; nrgn=N+1; 
save('myHCPcft.mat','nmrgn','nmrgn1','irgn','nrgn','-append');

%% 
sclr{1}=10; nmsclr1{1}='C-L'; 
sclr{2}=11; nmsclr1{2}='C-R'; 
sclr{3}=20; nmsclr1{3}='T-L'; 
sclr{4}=21; nmsclr1{4}='T-R';
sclr{5}=14; nmsclr1{5}='H-L'; 
sclr{6}=15; nmsclr1{6}='H-R';
sclr{7}=5; nmsclr1{7}='A-L'; 
sclr{8}=6; nmsclr1{8}='A-R';
sclr{9}=[7 12 16]; nmsclr1{9}='BD-L';  
sclr{10}=[7 13 17]; nmsclr1{10}='BD-R'; 
sclr{11}=[8 18 3]; nmsclr1{11}='S-L';
sclr{12}=[9 19 4]; nmsclr1{12}='S-R';
sclr=sclr'; nmsclr1=nmsclr1';

ep1=e.pos(ivx,:); 
N=length(sclr); ivxsclr=cell(N,1); pvxsclr=cell(N,1); 
for i=1:N
    ind=[]; 
    for j=1:length(sclr{i})
        s=sclr{i}(j); I=find(es1==s);
        if s==7, x=ep1(:,1); 
            if i==9, I(x(I)>=0)=[]; else, I(x(I)<0)=[]; end
        end
        ind=[ind; single(I)]; 
    end
    ivxsclr{i}=ind; pvxsclr{i}=single(ep1(ind,:));
end
save('myHCPcft.mat','sclr','nmsclr1','-append')

ivxclr={single(find(es1==1)); single(find(es1==2))};
irgnlr=[ivxclr; ivxsclr]; nmrgnlr1=[{'Cx-L'}; {'Cx-R'}; nmsclr1]; 
save('myHCPcft.mat','irgnlr','nmrgnlr1','-append');
