clear all
clc

cm=c1(:,1:3);
zm=c1(:,4);

xy=O(1:end,1:2)/10^5;
L=max(pdist(xy));
H=max(cm(:,3));
D=(L^2+H^2)^0.5; %max distance
m=7;
cl=linspace(0,L*0.8,m)'; % �ھڸ�ƪ��̤j����Z���ﱼ��ӭ�   
cz=linspace(0,30,m)';
      
% ����nan��
nanidx=find(isnan(zm));
zm(nanidx)=[];
cm(nanidx,:)=[];
cm(:,1:2)=cm(:,1:2)/10^5;
%��zzm 1-3�d  4-8��  9-13�t��
zm(zm==1|zm==2|zm==3)=1;
zm(zm==4|zm==5|zm==6|zm==7|zm==8)=2;
zm(zm==9|zm==10|zm==11|zm==12|zm==13)=3;
% remove 760m �[�`��   (ID WK-1E)
id760=find(cm(:,3)==760);
cm(9369:id760,:)=[];
zm(9369:id760,:)=[];
% O(33,:)=[];

%original data size = 9562
co=cm; 
zo=zm;

%% XY covariance

        method='kron';
        [d,P,o]=probatablecalc(cm(:,1:2),zm(:),cl,method);
        figure;
        probatableplot(d,P);
        

dfit=linspace(0,max(d)*0.9,m)'; % �令cl�̤j�Ȫ�90%
kstd=0.5*max(d)/length(d);
options=[1 2 1];
figure;
[PfitXY]=TPprobatablefit(dfit,d,P,o,kstd,options);

%% z covariance  �U���ӧO�� �A����
        method='kron';
for i=1:49  % �`�@�̦h49�f��
        idx=find(cm(:,1)==O(i,1)/10^5);
        [d(:,i),P,o(:,i)]=probatablecalc(cm(idx,3),zm(idx),cz,method);
        Ptemp11(:,i)=P{1,1};
        Ptemp12(:,i)=P{1,2};
        Ptemp13(:,i)=P{1,3};
        Ptemp22(:,i)=P{2,2};
        Ptemp23(:,i)=P{2,3};
        Ptemp33(:,i)=P{3,3};
        
end

        AVGd=nanmean(d,2);
        AVGo=nanmean(o,2);
        AVGP{1,1}=nanmean(Ptemp11,2);
        AVGP{1,2}=nanmean(Ptemp12,2);
        AVGP{2,1}=nanmean(Ptemp12,2);
        AVGP{1,3}=nanmean(Ptemp13,2);
        AVGP{3,1}=nanmean(Ptemp13,2);
        AVGP{2,2}=nanmean(Ptemp22,2);
        AVGP{2,3}=nanmean(Ptemp23,2);
        AVGP{3,2}=nanmean(Ptemp23,2);
        AVGP{3,3}=nanmean(Ptemp33,2);
         
        figure;
        probatableplot(AVGd,AVGP);
        
        dfit2=linspace(0,max(AVGd)*0.9,m)'; % �令cl�̤j�Ȫ�90%
        kstd=0.5*max(AVGd)/length(AVGd);
        options=[1 2 1];
        figure;
        [PfitZ]=TPprobatablefit(dfit2,AVGd,AVGP,AVGo,kstd,options);


%% �Ԧ�covariance

% z:XY=30:0.15
% z/200

cm(:,3)=cm(:,3)/200;

%%
       method='kron';
       [d,P,o]=probatablecalc(cm(:,1:3),zm(:),cl,method);
       figure;
       probatableplot(d,P);
        
dfit=linspace(0,max(d)*0.9,m)'; % �令cl�̤j�Ȫ�90%
kstd=0.5*max(d)/length(d);
options=[1 2 1];
figure;
[Pfit]=TPprobatablefit(dfit,d,P,o,kstd,options);

%% VVVVVVVV
crossO=[];
crossO=num2cell(O);

for i=1:49  % �`�@�̦h49�f��
    clear idn cn zn   
    idn=find(cm(:,1)==O(i,1)/10^5);
        
%% remove 
ccv=[];
zcv=[];
ccv=cm(idn,:);
zcv=zm(idn,:);
cn=cm;
zn=zm;
cn(idn,:)=[];
zn(idn,:)=[];

%%
nmax=10; % ��X�ӾF���I�ӱ���...�ݸ���I�K�צөw
dmax=10; % �̻��Ω�������F���I�Z��...�ݸ���I�K�צөw
options=[1 1e-1]; % �P�D�o����צ����A�V��V�C
[pkh]=BMEcatHard(ccv,cn(:,1:3),zn,dfit,Pfit,nmax,dmax,options);

%%
[trash,zkhmax]=max(pkh');
Result=[];
Result=[ccv zkhmax'];
crossvalidation=[];
crossvalidation=[Result zcv];

Matrix=[];
for k=1:3 %�������G
    clear m4 M
    m4=find(crossvalidation(:,4)==k);
    M=crossvalidation(m4,5);
    for l=1:3 %��l���
        Matrix(k,l)=size(find(M(:,1)==l),1);
    end
end

AXAX=crossvalidation(:,4)-crossvalidation(:,5);

crossO(i,4)=mat2cell(crossvalidation);
crossO(i,5)=mat2cell(Matrix);
crossO(i,6)=num2cell(size(find(AXAX(:,1)==0),1));
crossO(i,7)=num2cell(size(ccv,1));
crossO(i,8)=num2cell(size(find(AXAX(:,1)==0),1)/size(ccv,1));

end
%% VVVVVVVVVV
% Category01=Result((find(Result(:,4)==1)),:);
% Category02=Result((find(Result(:,4)==2)),:);
% Category03=Result((find(Result(:,4)==3)),:);

% scatter3(Category01(:,1),Category01(:,2),Category01(:,3),'rs');
% hold on;
% scatter3(Category02(:,1),Category02(:,2),Category02(:,3),'bs');
% scatter3(Category03(:,1),Category03(:,2),Category03(:,3),'ms');

save TP1m_xy&z_cross_all
