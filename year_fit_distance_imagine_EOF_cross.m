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

%整理harddata為ps各類別機率
for i=1:size(harddata,1)
    
    for j=1:3
        if  harddata(i,4)==j
            harddata(i,4+j)=1;
        else 
            harddata(i,4+j)=0;
        end
    end
end

harddata(:,3)=harddata(:,3)/200;%harddata是原始坐標需拉伸,softdata的坐標本身已拉伸
imaginedata(:,3)=imaginedata(:,3)/200;%imaginedata是原始坐標深度需拉伸
EOFdata(:,3)=EOFdata(:,3)/200;

cs=[harddata(:,1:3);softdata(:,1:3)];

ps=[harddata(:,5:7);softdata(:,6:8)];

%% VVVVVVVV
O(33,:)=[];
crossO=[];
crossO=num2cell(O);

for i=36:49 %:49  % 總共最多49口井
    clear idn cn zn   
    idn=find(cs(:,1)==O(i,1)/10^5);
        
%% remove 先抓出hard data要推估的點，再把虛擬井加進去
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


nmax=10; % 選幾個鄰近點來推估...看資料點密度而定
dmax=10; % 最遠用於推估的鄰近點距離...看資料點密度而定
options=[1 1e-2]; % 與求得的精度有關，越精越慢
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
for k=1:3 %推估結果
    clear m4 M
    m4=find(crossvalidation(:,4)==k);
    M=crossvalidation(m4,5);
    for l=1:3 %原始資料
        Matrix(k,l)=size(find(M(:,1)==l),1);
    end
end

AXAX=crossvalidation(:,4)-crossvalidation(:,5);

crossO(i,4)=mat2cell(crossvalidation); %每個點的推估結果
crossO(i,5)=mat2cell(Matrix); %準確率的9x9
crossO(i,6)=num2cell(size(find(AXAX(:,1)==0),1)); %類別推估正確的數量
crossO(i,7)=num2cell(size(ccv,1)); %該井的全部資料量
crossO(i,8)=num2cell(size(find(AXAX(:,1)==0),1)/size(ccv,1)); %準確率


save year_fit_distance_imagine_EOF_cross2

% %% 畫透視剖面對照圖（自己選定一個切面）
% Y=find(Result(:,2)==27.70); %選定北緯27.70
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

% 繪圖Betarange
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


