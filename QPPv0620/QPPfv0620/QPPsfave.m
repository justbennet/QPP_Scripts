function [QPP,TMX,nmx,CWQPP,C,QPPg2s,Tg0,isref,cij,tshij,sgnij]=...
    QPPsfave(QPPsa,B,tsh,PLc,nitrg,cth,PLh,G2Y)
[nsbj,nscn]=size(B); [nX,nt]=size(B{1,1}); PL=length(PLc); nXL=nX*PL; 
nT=nt*nscn; PLe=PL+sum(PLh); 

%% Comparing QPPas and finding ref QPPa (most similar to others)
Tj=zeros(nXL,nsbj,'single'); % To store all flattened normalized QPPas
Ti=cell(nsbj,1); Ti(:)={zeros(2*tsh+1,nXL,'single')}; % ~ all shifted QPPas
for is=1:nsbj
    T=QPPsa{is}; 
    if any(T), Tj(:,is)=Tcfn(T,PLc); Ti{is}=Tshcfn(T,tsh,PLc); end
end; clear T

z=zeros(nsbj,'single'); cij=z; tshij=z; sgnij=z; clear z
for is=1:nsbj
    c0=Ti{is}*Tj; 
    [cij(is,:),tshij(is,:)]=max(abs(c0));    
    for j=1:nsbj, sgnij(is,j)=sign(c0(tshij(is,j),j)); end
end; clear c0 Ti Tj
    
[~,isref]=max(sum(cij)); 

%% Timeshifting QPPas to ref QPPa & averaging
Tg0=zeros(nX,PLe,'single'); numT=0;
for is=1:nsbj
    T=QPPsa{is}; 
    if any(T)
        s=tshij(is,isref);
        if s<tsh+1, s=tsh+1-s; T=[T(:,s+1:end) zeros(nX,s,'single')]; 
        elseif s>tsh+1, s=s-(tsh+1); T=[zeros(nX,s,'single') T(:,1:end-s)];
        end
        Tg0=Tg0+T*sgnij(is,isref); numT=numT+1;
    end
end; clear QPPsa
Tg0=Tg0/numT; 

%% Correlating ave of QPPas with scans, finding maxima & obtaining GrpQPP
Tg=Tg0;
for itr=1:nitrg
    P=Tcfn(Tg,PLc);
    Tg=zeros(nX,PLe,'single'); TMX=cell(nsbj,1); nmx=zeros(nsbj,1); 
    C=zeros(nsbj,nT,'single'); QPPg2s=cell(nsbj,1);
for is=1:nsbj
    D=zeros(nX,nT,'single'); 
    for i=1:nscn, D(:,(i-1)*nt+1:i*nt)=B{is,i}; end
    warning off; [C(is,:),tmx]=Tcorr(D,[],P,nscn,1,cth); warning on;
    n=length(tmx); TMX{is}=tmx; nmx(is)=n;  
    if n, Ts=Tbld(D,tmx,PL,PLh,0); Tg=Tg+Ts; QPPg2s{is}=Ts/n; end
end; Tg=Tg/sum(nmx);
end; QPP=Tg; clear P Ts D

%% Correlation within GrpQPP between Yeo's RSNs
CWQPP=nan;
if nargin>7
    nY=max(G2Y); T=zeros(nY,PL,'single'); 
    for i=1:nY, T(i,:)=mean(QPP(G2Y==i,PLc),1); end; CWQPP=corr(T');
end
