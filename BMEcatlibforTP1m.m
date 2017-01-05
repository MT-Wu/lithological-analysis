%clear all
%clc

load c
load cm
load zm

xy=O(1:end,1:2);

% figure;
% colorplot(c(:,1:2),zh);
% axis square;
% xlabel('Easting');
% ylabel('Northing');
% title('Observed categories')

L=max(pdist(xy));
H=max(cm(:,3));
D=(L^2+H^2)^0.5; %max distance
m=7;
cl=linspace(0,D,m)'; % 根據資料的最大成對距離改掉兩個值     
       

%         method='kron';
%         [d,P,o]=probatablecalc(cm,zm,cl,method);
%         figure;
%         probatableplot(d,P);

load z1m_all.mat
        
dfit=linspace(0,max(d)*0.9,m)'; % 改成cl最大值的90%
kstd=0.5*max(d)/length(d);
options=[1 2 1];
figure;
[Pfit]=TPprobatablefit(dfit,d,P,o,kstd,options);

save Pfit 

minc=[2.9*10^5 2.762*10^6 0]; % oringin point
dc=[1000 1000 10]; % 切的網格長寬
nc=[18 2 76]; % 整個底圖切完後，X軸Y軸的格數
[ck]=creategrid(minc,dc,nc);

nmax=10; % 選幾個鄰近點來推估...看資料點密度而定
dmax=10000; % 最遠用於推估的鄰近點距離...看資料點密度而定
options=[1 1e-3]; % 與求得的精度有關，越精越慢

nanidx=find(isnan(zm));
zm(nanidx)=[];
cm(nanidx,:)=[];

[pkh]=BMEcatHard(ck,cm,zm,dfit,Pfit,nmax,dmax,options);

[trash,zkhmax]=max(pkh');
[ck1,ck2,zkhmax]=col2mat(ck,zkhmax');

% figure;
% pcolor(ck1,ck2,zkhmax);
% colormap hot
% hold on;
% shading flat
% colorplot(c,zh);
% axis square;
% xlabel('Easting');
% ylabel('Northing');
% title('Maximum probability estimates (hard data)');
% axis square;
% 
% save BMECATLIBtutorial
