function ttl=PLTUMET(nT,nscn,ylm,cth,met)
xlim([1 nT]); nt=nT/nscn; xticks(nt:nt:nscn*nT);
ylim(ylm); yticks([ylm(1) 0 cth ylm(2)]); grid on; box on;
if nargin>4
if any(met), cfr=[100 10]; met(1:2)=round(met(1:2).*cfr)./cfr;
    ttl=['C@maxima:' num2str(met(1)) ...
    '  \Deltatmaxima(s):' num2str(met(2)) '  #maxima:' num2str(met(3))]; 
else, ttl=[]; 
end
end