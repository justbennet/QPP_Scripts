clear; clc;
% % load('B_GWCR_HCPR40.mat','B'); nscn=size(B,2); 
% % is=1; D1=[]; for i=1:nscn, D1=[D1 B{is,i}]; end; save('SA_D1.mat','D1','nscn');
load('SA_D1.mat','D1','nscn'); nt=size(D1,2)/nscn;
PL=30; cth=[0.2 0.3]; ncth1=1; nitr=15; ssg=1; esg=nt-PL+1;
ITP=[]; for k=1:nscn, ITP=[ITP; (k-1)*nt+(ssg:esg)']; end
ITP=ITP(randperm(length(ITP))); ITP=ITP(1); e=5;

[b1TMX,a1QPP,c1itr]=f1run_qpp_algorithm(D1,nscn,PL,ITP,cth,ncth1,nitr);
% [b2TMX,a2QPP,c2itr]=f2QPP0620(D1,nscn,PL,ITP,cth,ncth1+1,nitr);
[b2TMX,a2QPP,c2itr]=f2QPP0918(D1,nscn,PL,ITP,cth,ncth1,nitr);

a0=corr(a1QPP(:),a2QPP(:));
[b3,i1,i2]=intersect(b1TMX,b2TMX); t1=b1TMX; t1(i1)=[]; t2=b2TMX; t2(i2)=[];

t1e=[]; for k=1:length(t1), t1e=[t1e; (t1(k)-e:t1(k)+e)']; end
[~,~,ii2]=intersect(sort(t1e),t2); b41=t2; b41(ii2)=[]; 
t2e=[]; for k=1:length(t2), t2e=[t2e; (t2(k)-e:t2(k)+e)']; end
[~,~,ii1]=intersect(sort(t2e),t1); b42=t1; b42(ii1)=[]; 

clear k i1 i2 ii1 ii2 t1 t2 t1e t2e ...
    nscn nt PL cth ncth1 nitr ssg esg ITP e D1
