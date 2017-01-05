function [B]=Betarange_MT()

%% 計算Beta的range(每個P的上下界）
%  資料參數接續正在fit的sample
%  先隨機抽取幾個樣本，算出n組probatablecalc
%  再把n組同樣位置的取出來bootstrap

Beta={};
Beta11=[]; Beta12=[]; Beta13=[];
Beta21=[]; Beta22=[]; Beta23=[];
Beta31=[]; Beta32=[]; Beta33=[];
Beta={Beta11 Beta12 Beta13;Beta21 Beta22 Beta23;Beta31 Beta32 Beta33};
Alldata=[cm(:,1:3) zm(:)]; %資料總數量

for n=1:100 %抽樣100次
    
    Prand=[];%隨機抽樣的P
    Prandfit=[];%抽樣後照dfit距離fit出來的P
   
    max=9000  %改成要抽樣的資料總數量
    A=randsample(1:9562,max); %函數預設不重複抽取，先亂數抽出所有資料的位置列A
    B=sortrows(A');  %將亂數抽出的位置點轉正為B列
    Prand=Alldata(B,1:4);

       method='kron';
       [d,P,o]=probatablecalc(Prand(:,1:3),Prand(:,4),cl,method);
       %figure;
       %probatableplot(d,P);

       %dfit要固定->和樣本實驗值對應

       %kstd=0.5*max(d)/length(d);
       options=[1 2 1];
       %figure;
       [Prandfit]=TPprobatablefit(dfit,d,P,o,kstd,options);
       %把這組抽樣的結果接到每個Beta表上（之後方便計算）
       Beta11=[Beta11 Prandfit{1,1}];
       Beta12=[Beta12 Prandfit{1,2}];
       Beta13=[Beta13 Prandfit{1,3}];
       Beta21=[Beta21 Prandfit{2,1}];
       Beta22=[Beta22 Prandfit{2,2}];
       Beta23=[Beta23 Prandfit{2,3}];
       Beta31=[Beta31 Prandfit{3,1}];
       Beta32=[Beta32 Prandfit{3,2}];
       Beta33=[Beta33 Prandfit{3,3}];
       Beta={Beta11 Beta12 Beta13;Beta21 Beta22 Beta23;Beta31 Beta32 Beta33};
       fprintf('Betaround: %d\n',round(n)); %看跑到第幾個round
       save Betatable Beta;
end

%Beta={Beta11 Beta12 Beta13;Beta21 Beta22 Beta23;Beta31 Beta32 Beta33};

%  抽完n組probatablecalc，把n組同樣的位置取出來做bootstrap
Beta_upper={};
Beta_lower={};

for i=1:3 
    for j=1:3
        Ci_lower=[];
        Ci_upper=[];
        for k=1:size(dfit,1)  %看前面抽樣幾次，先把資料叫出來
            data=Beta{i,j}(k,:); 
            Ci=[];
            Ci = bootci(5000,@mean,data); %直接算每個的信賴區間
            
            Ci_lower=[Ci_lower ; ci(1,:)];
            Ci_upper=[Ci_upper ; ci(2,:)];
        end
        
        Beta_lower{i,j}=Ci_lower;
        Beta_upper{i,j}=Ci_upper;
        save Betatable Beta Beta_lower Beta_upper;
    end
end

%找Beta上下界和實驗值的差，並選擇差較大的
Betarange={};
for i=1:3 
    for j=1:3
        lowerbound=[];
        upperbound=[];
        bound=[];
        Maxbound=[];
        
        lowerbound=Pfit{i,j}-Beta_lower{i,j};
        upperbound=Beta_upper{i,j}-Pfit{i,j};
        bound=[lowerbound(:,1) upperbound(:,1)];
        bound=roundn(bound*10^6,0); %後面找max需為整數
        Maxbound=max(bound,[],2); %找到較大的誤差範圍
        
        Betarange{i,j}=Maxbound/(10^6); %除回去
        
        save Betatable Beta Beta_lower Beta_upper Betarange;
    end
end

end
