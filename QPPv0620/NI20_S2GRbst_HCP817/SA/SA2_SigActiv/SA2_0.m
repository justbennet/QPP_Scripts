
%% Significant activation/deactivation within QPP
%% Part0) Random selection of segments to build null patterns
%%
clear; clc; nrep=50; p2s='SA2.mat'; p2='../../'; p1='../';
p2p=dir([p2 'Params_*.mat']); p2p=[p2 p2p.name];
load(p2p,'p2O','d2SA','nsbj','nscn','nt','PL','nT'); 
p2O=[p2 p2O]; d2SA=[p1 d2SA(3:end)];
load(p2O,'NMX'); nscng=nsbj*nscn;

%%
s=sum(NMX,3); % sr=round(s/10^3*10)/10;
nmxn=[max(s(:)) min(s(:))]; 
a=single(1:nt*nscng)'; 
ix=[]; for i=1:nscng, ix=[ix;(i-1)*nt+(nt-PL+2:nt)']; end
a(ix)=[]; l=length(a); % lr=round(l/10^6*10)/10;

TN=cell(nrep,2); TN(:)={cell(nsbj,1)};
for ir=1:nrep
    b=a(randperm(l));
    d=abs(diff(b)); I=find(d<=PL); ix=I;
    while any(I), d(I)=nan; I=find(d<=PL); ix=[ix;I]; end; clear d I
    b(ix)=[];
    I{1}=sort(b(1:nmxn(1))); I{2}=sort(b(end-nmxn(2)+1:end));
    for is=1:nsbj 
    for i=1:2    
       J=I{i}; J(J<(is-1)*nT+1|J>is*nT)=[]; TN{ir,i}{is}=J-(is-1)*nT;
    end
    end
end
save([d2SA p2s],'TN','nmxn');
