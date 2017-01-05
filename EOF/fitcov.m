function fitcov()
F=1;
%% Ozone��covariance  % from �ʭ�code
filename='covmap';
%%%%%%
% ���X���
load well_leveldata.mat
n=size(wellno,1); %�X�f��
% de-trend

trend=nanmean(waterlevel);
detrendlevel=waterlevel-repmat(trend,length(waterlevel),1);

% detrendlevel=waterlevel; % no de-trend baseline 
%
lonx=num(:,5);
laty=num(:,6);
tt=[1:length(waterlevel)]'*ones(1,n);
tt=reshape(tt',length(waterlevel)*n,1);

%%%%%%
cMS=[lonx laty];% x y�y��
tME=1:length(waterlevel); % t���ɶ��ǦC length(waterlevel) = 2005/2/1~2013/12/31

ch=[repmat(cMS,length(waterlevel),1) tt]; 
zh=reshape(detrendlevel',length(waterlevel)*n,1); % data

[sdmax,sdmin,tdmax,tdmin]=maxmind(cMS,tME);
[Z,cMS,tME,nanratio]=valstv2stg(ch,zh);

% ns= input('How do you want to lags of number(space) ? ');
ns=16;% �X��pair�I
ns=ns-1;
y= (sdmax/ns).*[1:ns-1];
para.rLag=[0,y,sdmax];
para.rLagTol=[0, para.rLag(2:end)- para.rLag(1:end-1)];
para.rLag=(para.rLag)';
                    
% nt= input('How do you want to lags of number(time) ? ');
tdmax=70;
% nt=15;
% nt=nt-1;
%y= (tdmax/nt).*[1:nt-1];
para.tLag=0:7:tdmax;
%[0,y,tdmax];
para.tLagTol=0.5*ones(size(para.tLag));
                     
[Cr,npr]=stcov(Z,cMS,tME,Z,cMS,tME,para.rLag,para.rLagTol...
                  ,para.tLag,para.tLagTol);
  
  Property={'Linewidth','Color'};
  Value={2,'b'};              
disp('Calculate Covariance finish')

%% ��covariance�Ѽư�

%by yourself to do
covmodelS={'exponentialC','exponentialC'};%change value
covparamS={[4.4 25377.34],[2.0 10000.56 ]}; %change value  
%            [A     B],[   C    D]   
covmodelT={'gaussianC','gaussianC'};%change value %gaussianC
covparamT={[4.4 150.0],[2.0 120]};%change value
%            [A     E],[   C    F]   
[covmodel,covparam,c,d]=covmodelST(covmodelS,covparamS,covmodelT,covparamT);

covmodel={'exponentialC/exponentialC','gaussianC/gaussianC'};
covparam={[4.4 25377.34 150.0],[4.4 10000.56 120]};
%          [A      B      E ],[ C    D    F]   

%% draw ���ν�
figure;
subplot(2,1,1); hold on;
plot(para.rLag(:,1),Cr(:,1),'or','MarkerSize',6,'Linewidth',1);
xlim([0 2.5*10^4])
%modelplot(0:0.1:rLag(end),covmodelS,covparamS,Property,Value);
modelplot(0:0.5:para.rLag(end),covmodelS,covparamS,Property,Value);
xlabel('Spatial lag (meter)');
ylabel('covariance ');
legend_handle=legend(gca,'show','Empirical covariance','Covariance model','Location','North');
set(legend_handle, 'Box', 'off');
set(legend_handle, 'Color', 'none');
%title('Spatial Covariance');
subplot(2,1,2); hold on;
plot(para.tLag(1,:),Cr(1,1:end),'or','MarkerSize',6,'Linewidth',1);
%modelplot(0:0.1:tLag(end),covmodelT,covparamT,Property,Value);
modelplot(0:0.5:para.tLag(end),covmodelT,covparamT,Property,Value);
xlabel('Time lag (day)');
ylabel('covariance');
% Create legend
legend_handle=legend(gca,'show','Empirical covariance','Covariance model','Location','North');
set(legend_handle, 'Box', 'off');
set(legend_handle, 'Color', 'none');
% title(filename(1:end-4)) ;
%title('Temporal Covariance');
% dirrefig=pwd;
com_fig=strcat(filename ,'.png');
saveas(figure(1) ,com_fig,'png');
