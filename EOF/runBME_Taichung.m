clear all
clc

function runBME()

% first run covdata to reatrive best covariance model
% and type here

%%  BME
nhmax=15;
nsmax=4;
order=NaN;
option=BMEoptions;
% dmax=[999999999 10000 10000]; % [spatial_range time_range spatial_range/time_range];
dmax=[20000 2000 1500];
%%%%% COV MODEL
covmodel={'exponentialC/exponentialC','gaussianC/gaussianC'};
covparam={[4.4 25377.34 150.0],[4.4 10000.56 120]};


%% 取出資料(先在UTF8下面跑好存到well_leveldata.mat)
load well_leveldata.mat
n=size(wellno,1); %幾口井
% load EOFtoLithlogy_normEOF
% [num txt]=xlsread('地下水觀測井資訊.xlsx','地下水觀測井資訊_交大提供分層資訊');
% wellno=txt(2:end,2);
% waterlevel=[];
% for i = 1:length(wellno)
%     
%     data=xlsread(strcat(wellno{i,1},'.txt'));
%     waterlevel=[waterlevel data];
%     % 由左至右順序為原讀入之順序
%     
% end

cd('20151022result')  %改存到新資料夾
%有NaN的話，要刪除有NAN的資料日期  從2009/1/1開始算
% de-trend
trend=nanmean(waterlevel);
detrendlevel=waterlevel-repmat(trend,length(waterlevel),1);
% detrendlevel=waterlevel; % no de-trend baseline

lonx=num(:,5);
laty=num(:,6);
tt=[1:length(waterlevel)]'*ones(1,n);
tt=reshape(tt',length(waterlevel)*n,1);

%%%%%%
cMS=[lonx laty];% x y座標
tME=1:length(waterlevel); % t的時間序列 length(waterlevel) = 2009/1/1~2014/12/31

ch=[repmat(cMS,length(waterlevel),1) tt]; 
zh=reshape(detrendlevel',length(waterlevel)*n,1); % data

%% creat grid file
% lat=linspace(min(laty),max(laty),10); %10*10網格
% lon=linspace(min(lonx),max(lonx),10);
% 
% gridnum=[];
% for i = 1:10
%     for j=1:10
%         gridnum=[gridnum ; lon(i) lat(j)];
%     end
% end
% %gridnum=coredata_wellsite(:,1:2);
% %save grid50by50 gridnum
% 
% grid=[];
% for t=1:length(waterlevel)
% temp=[gridnum ones(length(gridnum),1)*t];
% grid=[grid;temp];
% end
% gridnum=grid;
% save well98_t2191 gridnum

%% 
load well98_t2191

cs=[]; % soft data 的座標不能和hard data一樣
% zh=data(:,4);
ck1=gridnum;
nl=[]; % 長度是soft data 個數
limi=[];  % 要改
probdens=[];  % 要改
 zh2=repmat(trend',length(waterlevel),1);
 for i=1:size(ck1,1)
    [zk(i,:),info(i,:)]=BMEprobaMoments(ck1(i,:),ch,cs,zh,2,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,option);
    
    if isnan(zk(i,:))==1
        disp('stop')
    end
    if mod(i,98)==1
        fprintf('BMEestimation: %d percent\n',i./size(ck1,1)*100);
    end
%% 推trend
%         [zk2(i,:),info2(i,:)]=BMEprobaMoments(ck1(i,:),ch,cs,zh2,2,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,option);
%     
%     if isnan(zk2(i,:))==1
%         disp('stop')
%     end
%     
    % fprintf('trend estimation-%d/%d\n',i,size(ck1,1));
 end

 %% 加回trend
% newzk=zk+zk2;
 %% 不加回trend
 newzk=zk;
 save BMEresult_wellsite newzk ck1
%%如果有分層：
%  newzk=zk;
%  if F==1   % 第一層
%     save BMEresultF1 newzk ck1
%  else  % F==2　　 % 第二層
%     save BMEresultF2 newzk ck1
%  end
%% 畫底圖code，記得把GIS的資料夾放入所上的路徑上
% grid
%%
% gridnum=xlsread('grid200by200.xlsx');
%load BMEresult
% load grid50by50 % gridnum
% point=xlsread('台中市村里形心座標及各方法推估結果.xlsx');
% point=point(:,2:3);
% xi=point(:,1);
% yi=point(:,2);
% XI=gridnum(:,1);
% YI=gridnum(:,2);
% BMEestimate=zk(:,1);
% trend=xlsread('Ozone PM25各站整理平均.xlsx','O3');
% trend=trend(11,3:26);
% estimate=[];
% for t=1:24
% savefigname=strcat('ozone_',num2str(t));
% %curr_r=reshape(BMEestimate(200*200*(t-1)+1:200*200*t),200,200);%轉變成grid
% XI=reshape(XI,200,200);%轉變成grid
% YI=reshape(YI,200,200);%轉變成grid
% curr_r=BMEestimate(214*(t-1)+1:214*t);
% curr_r=curr_r+trend(t); % 把之前移掉的trend加回去
% [XI,YI,ZI] = griddata(xi,yi,curr_r,XI,YI,'linear');
% % 外插
% % index=find(isnan(ZI));
% % [ZI2] = griddata(xi,yi,curr_r,XI(index),YI(index),'nearest');
% % ZI(index)=ZI2;
% estimate=[estimate curr_r];
% 
% dirGIS='./GIS/';
% dirOut='./Mapping/';
% taipeifile=[dirGIS '台中市村里圖.shp'];
% taipei=shaperead(taipeifile);
%   figure; hold on;
%   [trash1,hh] = contourf(XI,YI,ZI,'LineStyle','none');
%   %contourf(xi,yi,ratiovalue,'LineStyle','none');
%   %surf(xi,yi,ratiovalue,'LineStyle','none')
%   cmap = hot; 
%   cmap = cmap(end:-1:1,:); colormap(cmap);
%   %caxis([min(min(ratiovalue)) max(max(ratiovalue))]);
%   
%   caxis([0 120]); 
%   %title(currnum(t,:));
% 
%  colorbar; 
%    for i=1:length(taipei)
%     XX1{i}=taipei(i,1).X;
%     YY1{i}=taipei(i,1).Y;
%     h=plot(XX1{i},YY1{i},'k--');
%     set(h,'LineWidth',2);
%    end
% % title(sprintf('PM2.5 on %d,%d',dtnum,name2_c));
%   axis image;
%   axis([205000,230566,2665000,2680000]);
%   % 北2692157.8 東230566.7 西184675.6 南2642452.8  for 台中
% %   axis([303475.946443,307966.515677,2766506.06272,2770720.356998]);%for Daon
%   hold on;
% % cd('savefloder')
%   saveas(gcf,savefigname,'png');%png
% end
% % close;
% % cd ..
% xlswrite('Ozone estimate',estimate);

