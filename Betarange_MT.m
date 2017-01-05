function [B]=Betarange_MT()

%% �p��Beta��range(�C��P���W�U�ɡ^
%  ��ưѼƱ��򥿦bfit��sample
%  ���H������X�Ӽ˥��A��Xn��probatablecalc
%  �A��n�զP�˦�m�����X��bootstrap

Beta={};
Beta11=[]; Beta12=[]; Beta13=[];
Beta21=[]; Beta22=[]; Beta23=[];
Beta31=[]; Beta32=[]; Beta33=[];
Beta={Beta11 Beta12 Beta13;Beta21 Beta22 Beta23;Beta31 Beta32 Beta33};
Alldata=[cm(:,1:3) zm(:)]; %����`�ƶq

for n=1:100 %���100��
    
    Prand=[];%�H����˪�P
    Prandfit=[];%��˫��dfit�Z��fit�X�Ӫ�P
   
    max=9000  %�令�n��˪�����`�ƶq
    A=randsample(1:9562,max); %��ƹw�]�����Ʃ���A���üƩ�X�Ҧ���ƪ���m�CA
    B=sortrows(A');  %�N�üƩ�X����m�I�ॿ��B�C
    Prand=Alldata(B,1:4);

       method='kron';
       [d,P,o]=probatablecalc(Prand(:,1:3),Prand(:,4),cl,method);
       %figure;
       %probatableplot(d,P);

       %dfit�n�T�w->�M�˥�����ȹ���

       %kstd=0.5*max(d)/length(d);
       options=[1 2 1];
       %figure;
       [Prandfit]=TPprobatablefit(dfit,d,P,o,kstd,options);
       %��o�թ�˪����G����C��Beta��W�]�����K�p��^
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
       fprintf('Betaround: %d\n',round(n)); %�ݶ]��ĴX��round
       save Betatable Beta;
end

%Beta={Beta11 Beta12 Beta13;Beta21 Beta22 Beta23;Beta31 Beta32 Beta33};

%  �⧹n��probatablecalc�A��n�զP�˪���m���X�Ӱ�bootstrap
Beta_upper={};
Beta_lower={};

for i=1:3 
    for j=1:3
        Ci_lower=[];
        Ci_upper=[];
        for k=1:size(dfit,1)  %�ݫe����˴X���A�����ƥs�X��
            data=Beta{i,j}(k,:); 
            Ci=[];
            Ci = bootci(5000,@mean,data); %������C�Ӫ��H��϶�
            
            Ci_lower=[Ci_lower ; ci(1,:)];
            Ci_upper=[Ci_upper ; ci(2,:)];
        end
        
        Beta_lower{i,j}=Ci_lower;
        Beta_upper{i,j}=Ci_upper;
        save Betatable Beta Beta_lower Beta_upper;
    end
end

%��Beta�W�U�ɩM����Ȫ��t�A�ÿ�ܮt���j��
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
        bound=roundn(bound*10^6,0); %�᭱��max�ݬ����
        Maxbound=max(bound,[],2); %�����j���~�t�d��
        
        Betarange{i,j}=Maxbound/(10^6); %���^�h
        
        save Betatable Beta Beta_lower Beta_upper Betarange;
    end
end

end
