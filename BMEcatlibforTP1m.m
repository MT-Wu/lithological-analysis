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
cl=linspace(0,D,m)'; % �ھڸ�ƪ��̤j����Z���ﱼ��ӭ�     
       

%         method='kron';
%         [d,P,o]=probatablecalc(cm,zm,cl,method);
%         figure;
%         probatableplot(d,P);

load z1m_all.mat
        
dfit=linspace(0,max(d)*0.9,m)'; % �令cl�̤j�Ȫ�90%
kstd=0.5*max(d)/length(d);
options=[1 2 1];
figure;
[Pfit]=TPprobatablefit(dfit,d,P,o,kstd,options);

save Pfit 

minc=[2.9*10^5 2.762*10^6 0]; % oringin point
dc=[1000 1000 10]; % ����������e
nc=[18 2 76]; % ��ө��Ϥ�����AX�bY�b�����
[ck]=creategrid(minc,dc,nc);

nmax=10; % ��X�ӾF���I�ӱ���...�ݸ���I�K�צөw
dmax=10000; % �̻��Ω�������F���I�Z��...�ݸ���I�K�צөw
options=[1 1e-3]; % �P�D�o����צ����A�V��V�C

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
