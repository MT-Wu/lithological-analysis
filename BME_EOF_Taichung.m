function BME_EOF_Taichung()
varimax_rot=0; % EOF varimax rotation : 0=off, 1=on
zscore_open=1; % nomarlization : 0=off, 1=on

% 地調所計畫  台中地下水 12口第一層井 
% 2009-2014年 daily
% 取出各年分各站的水位資料，
% 先BME 並進行EOF
% 測站順序依照:淺層水井list.xlsx
% 交大已內插完成
% xx: 第幾個EOF
% %% extract data

% 取出資料
load well_leveldata.mat

%% perform BME
% fitcov()
% runBME()
load BMEresult.mat   % ck1 newzk
waterlevel=reshape(newzk(:,1),100,2191);
waterlevel=waterlevel';

%% perform EOF
disp('start EOF');
% data standardize
% PCA
switch zscore_open
    case 0   % without normaliztion
        [COEFF,SCORE,latent,tsquare] = princomp(waterlevel(:,:));
    case 1  % with normaliztion
        [COEFF,SCORE,latent,tsquare] = princomp(zscore(waterlevel(:,:)));
end
percent=latent./sum(latent)*100;

EC=[];
EOF=[];
Lambda=[];

for xx=1:5
%% 時間是EC
h1=figure;
plot(SCORE(:,xx)) % <--
EC1=strcat('EOF',num2str(xx),' :  ',num2str(percent(xx)),'%'); % <-- 

title(EC1,'FontSize',12,'FontWeight','bold') %自行抬頭調整顯示的大小
EC=[EC SCORE(:,xx)];
ylabel('normalized EOF')  %改圖的標題
xlabel('Time')  %改圖的時間  %days from 2009/1/1
set(gca,'XTick',linspace(1,2191,7)); %2191是日資料的長度
set(gca,'XTickLabel',num2cell(2009:2015))

%% EOF1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% varimax rotation %%%%%
if varimax_rot==1
    rotatedCOEFF = rotatefactors(COEFF(:,1:5),'Method','varimax');
    EOF1=rotatedCOEFF(:,xx);% EOF1  <--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    EOF1=COEFF(:,xx);% EOF1
end
EOF=[EOF EOF1];
% lonx=num(:,5);
% laty=num(:,6);
lonx=ck1(1:100,1);
laty=ck1(1:100,2);
% xi = linspace(min(lonx),max(lonx),100);
% yi = linspace(min(laty),max(laty),100);
xi = linspace(210931.541,220854.606,100);  %對應到資料點的邊界
yi = linspace(2638534.234,2684037.067,100);

COV=var(COEFF(:,xx));

[XI,YI] = meshgrid(xi,yi);
Z=griddata(lonx,laty,EOF1,XI,YI,'linear');  %'v4'
% Z=griddata(lonx,laty,zeros(size(EOF1,1),size(EOF1,2)),XI,YI,'v4');
h2=figure;
contourf(XI,YI,Z,100,'linestyle','none');
hold on;
colorbar;
% plot(lonx,laty,'rs','MarkerFaceColor','g','MarkerSize',10)
% text(lonx,laty,wellno(1:end,1),'FontWeight','bold')
plot(num(:,5),num(:,6),'rs','MarkerFaceColor','g','MarkerSize',10)
text(num(:,5),num(:,6),wellno(1:end,1),'FontWeight','bold','FontSize',12)
% caxis([-0.5 0.5]);
% colormap(flipud(hot));
colormap(jet);
title(EC1,'FontSize',12,'FontWeight','bold');
xlabel('longitude');
ylabel('latitude');
mapplot();
axis equal
axis tight

%axis([200000 230000 2638534.234 2684037.067]); %圖上點的邊界
axis([210931.541 220854.606 2638534.234 2684037.067]); %圖上點的邊界

figname1=strcat('EC',num2str(xx));
% saveas(h1,figname1,'fig');
saveas(h1,figname1,'png');
figname2=strcat('EOF',num2str(xx));
% saveas(h2,figname2,'fig');
saveas(h2,figname2,'png');


end
%% write xls EC
% lambda
Lambda=latent(:,1);% Lambda100*1

xlswrite('Taichung_EC',EC);
xlswrite('Taichung_EOF',EOF);
xlswrite('Taichung_Lambda',Lambda);

function mapplot()
    dirGIS='./GIS/';
    dirOut='./20150522result/';
    shpfile=[dirGIS 'TaichungGW2D_TWD97.shp']; % TaichungGW2D_TWD97 TWN_RIVER TWN_COUNTY_97_geo
    world=shaperead(shpfile);
    disp('continue')
    
    for i=1:length(world)
        XX1{i}=world(i,1).X;
        YY1{i}=world(i,1).Y;
        h=plot(XX1{i},YY1{i},'k');
        idx=find(XX1{i}<0);
        XX1{i}(idx)=XX1{i}(idx)+360;
        if ~(max(XX1{i}-min(XX1{i})>350))
            h=plot(XX1{i},YY1{i},'k');
        end
        set(h,'LineWidth',1);
    end
    
    