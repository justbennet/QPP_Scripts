function [dFC,hFC,pFC,qFCN,hFCN]=QPPsdFC(B,FCr,ibY,iG2Y,fcbn,fcth)
[nsbj,nscn]=size(B); [nX,nt]=size(B{1,1}); nY=length(iG2Y); nP=length(FCr);

%% FC of original scan
cz=zeros(nX);
for is=1:nsbj, czs=zeros(nX);
for iscn=1:nscn
    d=zeros(nX,nt);
    for i=1:nY, d(ibY(i)+1:ibY(i+1),:)=B{is,iscn}(iG2Y{i},:); end
    c=corr(d'); c(c==1)=0.9999; czs=czs+0.5*log((1+c)./(1-c));
end; cz=cz+czs/nscn; 
end; cz=cz/nsbj; FC=(exp(2*cz)-1)./(exp(2*cz)+1); FC(FC>=0.9999)=1; 

%% Combining FC of the original & residual scans as upper & lower triangles
dFC=cell(nP,1); dFC(:)={nan(nX,'single')}; 
iu=find(triu(ones(nX),1)); il=find(tril(ones(nX),-1));
for ip=1:nP, dFC{ip}(iu)=FC(iu); dFC{ip}(il)=FCr{ip}(il); end; clear iu il

%% Histogram of FC values (hFC) & percentage of FC values > fcth
FC(eye(nX)>0)=nan; FC=FC(:); 
hFC=zeros(1+nP,length(fcbn),'single');
hFC(1,:)=hist(FC,fcbn);
pFC=zeros(1+nP,3,'single'); nX2=nX^2-nX; 
pFC(1,:)=[length(find(FC<-fcth)) ...
    length(find(FC>fcth)) length(find(abs(FC)>fcth))]/nX2*100; clear FC
for ip=1:nP
    A=FCr{ip}; A(eye(nX)>0)=nan; A=A(:); 
    hFC(1+ip,:)=hist(A,fcbn);
    pFC(1+ip,:)=[length(find(A<-fcth)) ...
        length(find(A>fcth)) length(find(abs(A)>fcth))]/nX2*100; 
end; clear A

%% Null distribution for FC
cz=zeros(nX);
for is=1:nsbj, czs=zeros(nX);
for iscn=1:nscn
    amp=abs(fft(B{is,iscn},[],2)); 
    r=rand(nX,nt,'single'); r=zscore(r,[],2); pha=angle(fft(r,[],2));  
    N=ifft(amp.*exp(sqrt(-1)*pha),[],2); N=real(N);
    d=zeros(nX,nt);
    for i=1:nY, d(ibY(i)+1:ibY(i+1),:)=N(iG2Y{i},:); end
    c=corr(d'); c(c==1)=0.9999; czs=czs+0.5*log((1+c)./(1-c));
end; cz=cz+czs/nscn;
end; cz=cz/nsbj; FCN=(exp(2*cz)-1)./(exp(2*cz)+1); FCN(FCN>=0.9999)=1;
FCN(eye(nX)>0)=nan; FCN=FCN(:); hFCN=single(hist(FCN,fcbn));
qFCN=quantile(FCN,0.99);
