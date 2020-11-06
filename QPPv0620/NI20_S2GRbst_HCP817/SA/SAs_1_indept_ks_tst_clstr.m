
%% Testing significance of time differences btwn QPP's cluster per region
%% using independent t-test and ks-test, without matching group sizes
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=p2p.name;
load([p1 p2p],'p2V','nP','nSD','nsd','nclst','nrgn','nmrgn1','d2SAplt'); 
nth=15; nrc=nrgn*nclst; benf=nrc*(nrc-1)/2; alph=0.01/benf;
load([p1 p2V],'QPPv2clst_te'); 
ss=get(0,'Screensize'); nmte='pd';

%% Common part to t-/ks-test
nte=nan(nrgn,nclst,2,nP,nSD,'single');
for ip=1:nP
for isd=1:nSD
for ig=1:nrgn
for ic=1:nclst
    te=QPPv2clst_te{ip,isd}{ig,ic};
    if ~isempty(te)
        for ie=1:2
            t=te(:,ie); t(isnan(t))=[];
            nte(ig,ic,ie,ip,isd)=length(t);
        end
    end
end
end
end
end; clear te t
nte(nte<nth)=nan;

icR=cell(nrgn,2,nP,nSD);
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:2
for ig=1:nrgn
    icR{ig,ie,ip,isd}=find(~isnan(nte(ig,:,ie,ip,isd)));
end
end
end
end

%% uncomment t-/ks-test whichever is to test & plot
pttst_te=nan(nrc,nrc,2,nP,nSD,'single');
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:2
for ig1=1:nrgn
    icr1=icR{ig1,ie,ip,isd}; 
    if ~isempty(icr1)
    for ic1=icr1
        te1=QPPv2clst_te{ip,isd}{ig1,ic1}(:,ie); 
        te1(isnan(te1))=[];
        for ig2=1:nrgn
            icr2=icR{ig2,ie,ip,isd};
            if ~isempty(icr2)
            for ic2=icr2
                te2=QPPv2clst_te{ip,isd}{ig2,ic2}(:,ie);
                te2(isnan(te2))=[];
                [h,p]=ttest2(te1,te2,'Alpha',alph,'vartype','unequal');
%                 [h,p]=kstest2(te1,te2,'Alpha',alph);
                i=(ig2-1)*nclst+ic2; j=(ig1-1)*nclst+ic1;
                if ~h, pttst_te(i,j,ie,ip,isd)=1;
                else, pttst_te(i,j,ie,ip,isd)=p;
                end
            end
            end
        end
    end
    end
end
end
end
end

nrc=nrgn*nclst; il=nclst/2:nclst:nrc; l=cell(nrc,1); 
for i=1:length(il), l{il(i)}=nmrgn1{i}; end
for ip=1:nP
for isd=1%:nsd(ip)
for ie=1%:2
    figure; set(gcf,'Position',ss); 
    p=pttst_te(:,:,ie,ip,isd); 
    im=imagesc(p); set(im,'AlphaData',~isnan(p)); colorbar; axis square; 
    xticks(1:nrc); yticks(1:nrc); xticklabels(l); yticklabels(l); 
    hold on; grid on; 
    for i=nclst+0.5:nclst:nrc-nclst+0.5, plot([i i],[0 nrc+1],'k'); end
    for i=nclst+0.5:nclst:nrc-nclst+0.5, plot([0 nrc+1],[i i],'k'); end
    saveas(gcf,[d2SAplt 'SAs_1_indep-ttst_QPP' num2str(ip) ...
        'SD' num2str(isd) '_' num2str(ie) 't' nmte(ie) '.png']); close
    %     saveas(gcf,[d2SAplt 'SAs_1_ks-ttst_QPP' num2str(ip) ...
%     'SD' num2str(isd) '_' num2str(ie) 't' nmte(ie) '.png']); close
end
end
end

