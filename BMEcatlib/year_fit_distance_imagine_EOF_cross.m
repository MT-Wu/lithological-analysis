clear all
clc

%%
load zm;
load harddata;
load softdata;
load imagine_data;
load origin_data;
load distance_pdf;
load EOFdata;
load Betatable;

%��zharddata��ps�U���O���v
for i=1:size(harddata,1)
    
    for j=1:3
        if  harddata(i,4)==j
            harddata(i,4+j)=1;
        else 
            harddata(i,4+j)=0;
        end
    end
end

harddata(:,3)=harddata(:,3)/200;%harddata�O��l���лݩԦ�,softdata�����Х����w�Ԧ�
imaginedata(:,3)=imaginedata(:,3)/200;%imaginedata�O��l���в`�׻ݩԦ�
EOFdata(:,3)=EOFdata(:,3)/200;

cs=[harddata(:,1:3);softdata(:,1:3)];

ps=[harddata(:,5:7);softdata(:,6:8)];

%% VVVVVVVV
O(33,:)=[];
crossO=[];
crossO=num2cell(O);

for i=36:49 %:49  % �`�@�̦h49�f��
    clear idn cn zn   
    idn=find(cs(:,1)==O(i,1)/10^5);
        
%% remove ����Xhard data�n�������I�A�A��������[�i�h
ccv=[];
zcv=[];
ccv=cs(idn,:);
zcv=ps(idn,:);
cn=cs;
zn=ps;
cn(idn,:)=[];
zn(idn,:)=[];
cn=[cn(:,1:3);imaginedata(:,1:3);EOFdata(:,1:3)];
zn=[zn(:,1:3);imaginedata(:,5:7);EOFdata(:,11:13)];


nmax=10; % ��X�ӾF���I�ӱ���...�ݸ���I�K�צөw
dmax=10; % �̻��Ω�������F���I�Z��...�ݸ���I�K�צөw
options=[1 1e-2]; % �P�D�o����צ����A�V��V�C
[pkh]=BMEcatPdf_MT(ccv,cn(:,1:3),zn(:,1:3),dfit,Pfit,Betarange,nmax,dmax,options);
%  [pkh]=BMEcatPdf(ccv,cn(:,1:3),zn(:,1:3),dfit,Pfit,nmax,dmax,options);

%%
[trash,zkhmax]=max(pkh');
[trash2,zcvmax]=max(zcv');
Result=[];
Result=[ccv zkhmax'];
crossvalidation=[];
crossvalidation=[Result zcvmax'];

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

crossO(i,4)=mat2cell(crossvalidation); %�C���I���������G
crossO(i,5)=mat2cell(Matrix); %�ǽT�v��9x9
crossO(i,6)=num2cell(size(find(AXAX(:,1)==0),1)); %���O�������T���ƶq
crossO(i,7)=num2cell(size(ccv,1)); %�Ӥ���������ƶq
crossO(i,8)=num2cell(size(find(AXAX(:,1)==0),1)/size(ccv,1)); %�ǽT�v


save year_fit_distance_imagine_EOF_cross2

% %% �e�z���孱��ӹϡ]�ۤv��w�@�Ӥ����^
% Y=find(Result(:,2)==27.70); %��w�_�n27.70
% X=Result(Y,1);
% Z=Result(Y,3);
% C=Result(Y,4);
% C1=find(C==1);
% C2=find(C==2);
% C3=find(C==3);
% 
% scatter(X(C1),Z(C1),'rs');
% hold on;
% set(gcf,'color','none');
% set(gca,'color','none'); 
% scatter(X(C2),Z(C2),'bs');
% scatter(X(C3),Z(C3),'ms');

% ø��Betarange
% figure;
% probatableplot(dfit,Pfit);
% hold on;
% probatableplot(dfit,Beta_lower);
% 
% for i=1:3
%     for j=i+1:3
%         
x=Pfit{i,j};
y=Beta_lower{i,j};
z=Beta_upper{i,j};
% plot(dfit,x);
% hold on;
% plot(dfit,y,'m-');
% plot(dfit,z,'m-');
Range=[z x y];
% 
%     end
% end

end


