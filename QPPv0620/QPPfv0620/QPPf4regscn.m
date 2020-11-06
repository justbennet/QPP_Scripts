function [Dr,Cr,FCr]=QPPf4regscn(D,TT,CT,nscn,PL,PLc,ibY,iG2Y,nz)
[nX,nT]=size(D); nt=nT/nscn; nXL=nX*PL; npr=length(TT);
if size(TT{1},2)>PL, for i=1:npr, TT{i}=TT{i}(:,PLc); end; end

%%
Dr=zeros(nX,nT,'single'); Cr=zeros(npr,nT,'single'); FCr=zeros(nX);
for iscn=1:nscn
    %% Regressing QPPs 1 to npr
    tscn=(1:nt)+(iscn-1)*nt;
    t=tscn(PL:end); % << starting from PL due to valid option of conv
    for ix=1:nX
        r=zeros(nt-PL+1,npr,'single'); % << -PL+1 due to valid option
        for i=1:npr, r(:,i)=conv(CT(i,tscn),TT{i}(ix,:),'valid')'; end
        beta=(r'*r)\r'*D(ix,t)'; 
        Dr(ix,t)=D(ix,t)-(r*beta)';
    end; Dr(:,t)=zscore(Dr(:,t),[],2);

    %% Correlating QPPs 1-npr with residual scan (goodness of regression)
    for ipr=1:npr
        P=TT{ipr}(:); P=P-sum(P)/nXL; P=P'/sqrt(P'*P); 
    for isg=PL:nt-PL+1 % << starting from PL ~
        tsg=isg+(iscn-1)*nt; S=Dr(:,tsg:tsg+PL-1); 
        S=S(:); S=S-sum(S)/nXL; S=S/sqrt(S'*S); Cr(ipr,tsg)=P*S;
    end
    end   

    %% FC of residuals (FCr), averaged across scans using Fisher transform
    if nargin>6
        d=zeros(nX,nt-PL+1); % << -PL+1 ~
        for i=1:length(iG2Y), d(ibY(i)+1:ibY(i+1),:)=Dr(iG2Y{i},t); end
        c=corr(d'); c(c==1)=0.9999; FCr=FCr+0.5*log((1+c)./(1-c));
    end
end
if nargin>6, FCr=FCr/nscn;
    if nz, FCr=(exp(2*FCr)-1)./(exp(2*FCr)+1); FCr(FCr>=0.9999)=1; end
end

