function [c,sh]=Tcomp(T1,T2,PLc,tsh)
P1=Tshcfn(T1,tsh,PLc);
P2=Tcfn(T2,PLc);
c0=P1*P2';
[c,i]=max(abs(c0));
c=c*sign(c0(i));
sh=-tsh:tsh; sh=sh(i);