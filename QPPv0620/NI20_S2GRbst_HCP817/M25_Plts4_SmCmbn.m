
%% Combining screenshots of summary maps of QPPs (captured by workbench)
%%
clear; clc; close all;
nP=3; pe='.png'; ss=get(0,'screensize');

%% for the main figures
%%
ps='./Plts/cfts2/'; m={'clstr_','tp_','prcl'}; 
v={'_s_','_c_','_a_'}; nv=[2 2 1]; 
o={[1,7],[3,2];[8,0],[6,4];[5,0],[]};
v0='_C_L';

ip=1; im=1; iv=1; i=1;
f=[ps m{im} num2str(ip) v{iv} num2str(i) pe]; 
I=imread(f); % figure; y0=size(I,1); image(I); axis image
[y,x,~]=size(I); x2=x/2;

f=[ps m{im} num2str(ip) v0 pe]; 
I=imread(f); % figure; set(gcf,'Position',ss); image(I); axis image
y0=size(I,1); yp=25;
I=I(yp+1:end-yp,:,:); % figure; set(gcf,'Position',ss); image(I); axis image
I=imresize(I,[y nan]); % figure; set(gcf,'Position',ss); image(I); axis image

%% 
for ip=1:nP
for im=1:3 
    ps2=[ps m{im} num2str(ip)]; if im==3, ps2=[ps m{im}]; end
    J=zeros(y,x2*8,3,'uint8');
    for iv=1:3
    for in=1:nv(iv)
        f=[ps2 v{iv} num2str(in) pe]; I=imread(f); 
        for io=1:2
            j=o{iv,in}(io);
            if j, J(:,(j-1)*x2+1:j*x2,:)=I(:,(io-1)*x2+1:io*x2,:); end
        end
    end
    end
    f=[ps2 v0 pe];
    I=imread(f); I=I(yp+1:end-yp,:,:); I=imresize(I,[y nan]); 
    J=[I J]; imwrite(J,[ps2 pe],'Mode','lossless');
end
end

%% with borders for supplementary figures
%%
% ps='./Plts/cfts2/id/'; m={'','_tp','_td'}; n={'Y7_','G360_'};
% 
% in=1; ip=1; im=1; 
% f=[ps n{in} num2str(ip) m{im} pe]; 
% I=imread(f); % figure; image(I); axis image
% [y,x,~]=size(I);
% d=40;
% Zx=zeros(y*4+d,d,3,'uint8'); 
% Zy=zeros(d,x,3,'uint8');
% 
% %%
% for ip=1:nP
%     K=[];
%     for in=1:2
%         f=[ps n{in}(1:end-1) pe]; J=imread(f);
%         J=[J; Zy];
%         for im=1:3
%             f=[ps n{in} num2str(ip) m{im} pe]; I=imread(f); 
%             J=[J;I];
%         end
%         K=[K Zx J];
%     end
%     % figure; image(K); axis image
%     imwrite(K,[ps num2str(ip) pe],'Mode','lossless');
% end

