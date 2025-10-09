clc; close all;clear 
%% 파일 불러오기
load Temperature.mat % 수온 (섭씨)
% load Salt.mat % 염분 (psu)
% load Oxygen.mat % 용존산소량 (mL/L)
%% 트릭
A = T; % 수온
clearvars T
% 
% A = S; % 염분
% clearvars S
% 
% A = DO; % 용존산소량
% clearvars DO

%% set line number
B = A(A(:,1)==207,:);
clearvars A
%% 변수 저장
cruiseline = B(:,1);
station = B(:,2);
depth = B(:,3);
year = B(:,4);
month = B(:,5);
day = B(:,6);
data = B(:,7); % 수온, 염분, 용존산소량 등
data_qc = B(:,8);
%% 정점, 깊이, 행렬의 크기
station_num = unique(station); % 
depth_lev = unique(depth);
vert_temp = nan(length(depth_lev),480);
vert208 = nan(length(depth_lev),480*length(station_num));
depth208 = nan(length(depth_lev),480*length(station_num));
%% QC Flag
qc = data_qc==1|data_qc==2;
num_qc1 = zeros(40,12,length(station_num));
num_qc2 = zeros(40,12,length(station_num));
%% count QC
for i=1:40
    for j=1:12
        for k=1:length(station_num)
             if data_qc(year==(1979+i)&month==(j)&qc)
                num_qc1(i,j,k)=nnz(data_qc(year==(1979+i)&month==j&station==k&data_qc==1));
                num_qc2(i,j,k)=nnz(data_qc(year==(1979+i)&month==j&station==k&data_qc==2));
             else
                 continue
             end
        end
    end
end
for k=1:length(station_num)
    for i=1:40
        figure;
        histogram(num_qc1(i,:,k))
    end
end
%% 자료 저장
for n=1:length(station_num)
    for i=1:length(depth_lev) % 185
        for j=1:12 % 12개월
            for k=1:40 % 1980-2019년 
                if ~isnan(data(cruiseline==208&station==station_num(n)&year==(1979+k)&month==j&depth==depth_lev(i)&qc))
                        vert208(i,(j+12*(k-1)+480*(n-1)))=unique(data(cruiseline==208&station==station_num(n)&year==(1979+k)...
                            &month==j&depth==depth_lev(i)&qc));
                        depth208(i,(j+12*(k-1)+480*(n-1))) = unique(depth(cruiseline==208&station==station_num(n)&year==(1979+k)...
                            &month==j&depth==depth_lev(i)&qc));

                else
                    continue
                end
             end
        end
    end
end
numberofvert208 = nnz(~isnan(vert208));
%% array reshape
vert = reshape(vert208,length(depth_lev),[],length(station_num));
% v208_01 = squeeze(vert(:,2:end,1));
% v208_02 = squeeze(vert(:,2:end,2));
% v208_03 = squeeze(vert(:,2:end,3));
% v208_04 = squeeze(vert(:,2:end,4));
% v208_05 = squeeze(vert(:,2:end,5));
% v208_06 = squeeze(vert(:,2:end,6));
% v208_07 = squeeze(vert(:,2:end,7));
% v208_08 = squeeze(vert(:,2:end,8));
% v208_09 = squeeze(vert(:,2:end,9));
d208 = nan(1,127);
for i=1:127
    d208_temp = unique(depth208(i,:));
    d208(i) = d208_temp(~isnan(d208_temp));
end
depth208_re = reshape(depth208,length(depth_lev),[],length(station_num));
% d208_01 = squeeze(depth208_re(:,2:end,1));
%% interpolation
depth_max = [70 130 150 143 140 135 130 130 130];
depth_standard = [0 10 20 30 50 75 100 125 150 200 250 300 400 500]; % 14 표준 수층
v208_interp=[];
%% Depth interpolation
depth1 = 0:1:205; % maximum depth 150 m
for j=1:length(station_num)
    v208_i = squeeze(vert(:,2:end,j)); %
    d208_i = squeeze(depth208_re(:,2:end,j));
    v208_i_interp1 = nan(length(depth1),length(v208_i));
    depth_len = length(depth1);
    for i =1:size(v208_i,2)
        dmin = ceil(min(d208_i(:,i)));
        dmax = floor(max(d208_i(:,i)));
        depth_interp = dmin:1:dmax;
        i_nan = find(~isnan(v208_i(:,i)));
        i_nan2 = i_nan(~isnan(d208_i(i_nan,i)));
        if ~isnan(i_nan2)
            F = griddedInterpolant(d208_i(i_nan2,i),v208_i(i_nan2,i),'spline');
            v208_i_interp1(find(depth1==depth_interp(1)):find(depth1==depth_interp(end)),i)= ...
                F(depth_interp);
        else
            continue
        end
    end
    depth_standard_index = depth_standard <=depth_max(j);
    depth_level_index = depth_standard(depth_standard_index)+1;
    v208_i_interp1_mod =v208_i_interp1(depth_level_index,1:2:end);
%% Time interpolation
    v208_i_interp = nan(length(depth_level_index),492);
   for i = 1:size(v208_i_interp1_mod,1)
        i_nan = find(~isnan(v208_i_interp1_mod(i,:)));

            F = griddedInterpolant(i_nan,v208_i_interp1_mod(i,i_nan),"spline");
            v208_i_interp(i,1:492) = F(1:492);
   end    

% v208_01_interp = squeeze(v208_01_interp(:,1:2:end));

%% interpolation check
% After Depth interpolation
    for i=1:length(depth_level_index)
        figure;
        plot(v208_i(depth_lev==depth_standard(i),1:2:480),'or',LineWidth=1.0)
        hold on
        plot(v208_i_interp1_mod(i,1:492),'b-',LineWidth=1.0)
        hold off
        ylim([-10 40])
        xlim([1 241])
        title(strcat("208 Line","\_","0",num2str(j),"\_",num2str(depth_standard(i))," m"),"FontSize",16,"FontWeight","bold");
        xticks(1:60:241)
        xticklabels(1980:10:2020)
        legend("raw data","interpolated data")
        xlabel("time","FontSize",14,"FontWeight","bold")
        ylabel("Temperature ( \circC)","FontSize",14,"FontWeight","bold")
        print(strcat("Di207","_","0",num2str(j),"_",num2str(depth_standard(i)),"m"),"-depsc","-tiff")
        close(gcf)
    end
%     legend("raw data","interpolated data","Location","bestoutside")
%     title(tl,"Depth interpolation","FontSize",20,"FontWeight","bold")
%     xlabel(tl,"time","FontSize",14,"FontWeight","bold")
%     ylabel(tl,"Temperature ( \circC)","FontSize",14,"FontWeight","bold")

% After Time interpolation
    for i=1:length(depth_level_index)
        figure;
        plot(v208_i(depth_lev==depth_standard(i),1:2:480),'or',LineWidth=1.0)
        hold on
        plot(v208_i_interp(i,1:492),'b-',LineWidth=1.0)
        hold off
        ylim([-10 40])
        xlim([1 241])
        title(strcat("208 Line","\_","0",num2str(j),"\_",num2str(depth_standard(i))," m"),"FontSize",16,"FontWeight","bold");
        xticks(1:60:241)
        xticklabels(1980:10:2020)
        legend("raw data","interpolated data")
        xlabel("time","FontSize",14,"FontWeight","bold")
        ylabel("Temperature ( \circC)","FontSize",14,"FontWeight","bold")
        print(strcat("Ti207","_","0",num2str(j),"_",num2str(depth_standard(i)),"m"),"-depsc","-tiff")
        close(gcf)
    end
%     legend("raw data","interpolated data","Location","bestoutside")
%     title(tl,"Time interpolation","FontSize",20,"FontWeight","bold")
%     xlabel(tl,"time","FontSize",14,"FontWeight","bold")
%     ylabel(tl,"Temperature ( \circC)","FontSize",14,"FontWeight","bold")
    v208_interp = vertcat(v208_interp,v208_i_interp);
    save("line207","v208_interp")
end
whos
% clearvars v208_i v208_i_interp v208_i_interp1_mod vert vert208 vert temp depth_level_index depth_standard_index depth208 depth208_re d208_temp d208 d208_i vert_temp v208_i_interp1
%% calculate distance between the stations
lats = [35.475 35.405 35.2967 35.185 35.08 34.97 34.8667 34.755 34.6483]; % 북위
lons = [129.455 129.5633 129.715 129.8783 130.033 130.1983 130.3483 130.51 130.675]; % 동경
d = nan(1,9);
d(1) = 0;
for i=1:8
    d_temp = sqrt(lons(i+1)^2 + lats(i+1)^2) - sqrt(lons(i)^2 +lats(i)^2);
    d(i+1) = d_temp;
end
%% plot vertical section mean
mean208=squeeze(mean(vert208_re,2,'omitnan'));

figure(1)
contourf(cumsum(d),1:depth_len,mean208,-5:0.1:30,'LineStyle','none') % 수온 평균
% contourf(cumsum(d),1:depth_len,mean102,'LineStyle','none') % 염분평균
% contourf(cumsum(d),1:depth_len,mean102,'LineStyle','none') % 용존산소평균

set(gca,'YDir','reverse')
yticks(linspace(1,depth_len,6))
yticklabels(0:100:500)
xticks(cumsum(d))
xticklabels(station_num)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
xlabel("Station No.","FontSize",14,"FontWeight","bold")
grid on
colormap("jet")
colorbar
clim([0 25])
caxis([0 25])
title('102 정선 수온 평균 ( \circC) 1980 - 2019','FontSize',20,'FontWeight','bold')
% print('mean102T','-depsc','-tiff')

% colormap("parula")
% colorbar
% clim([32 34.5])
% caxis([32 34.5])
% title('102 정선 염분 평균 (psu) 1980 - 2019','FontSize',20,'FontWeight','bold')
% print('mean102S','-depsc','-tiff')

% colormap("parula")
% colorbar
% clim([0 10])
% caxis([0 10])
% title('102 정선 DO 평균 ( mL/L) 1980 - 2019','FontSize',20,'FontWeight','bold')
% print('mean102DO','-depsc','-tiff')
%% plot vertical section std
std102=squeeze(std(vert208_re,0,2,'omitnan'));
figure(2)
contourf(cumsum(d),1:depth_len,std102,0:0.1:8,'LineStyle','none') % 수온 표준편차 
% contourf(cumsum(d),1:depth_len,std102,'LineStyle','none') %염분 표준편차
% contourf(cumsum(d),1:depth_len,std102,'LineStyle','none') %용존산소량 표준편차

set(gca,'YDir','reverse')
yticks(linspace(1,depth_len,6))
yticklabels(0:100:500)
xticks(cumsum(d))
xticklabels(station_num)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
xlabel("Station No.","FontSize",14,"FontWeight","bold")
grid on
colormap("jet")
colorbar('Ticks',[0 2 4 6 8],'TickLabels',{'0','2','4','6','8'})
clim([0 8])
caxis([0 8])
title("102 정선 수온 표준편차 ( \circC) 1980 - 2019",'FontSize',20,'FontWeight','bold')
% print('std102T','-depsc','-tiff')
% 
% colormap("parula")
% colorbar('Ticks',[0 0.5 1.0 1.5 2],'TickLabels',{'0','0.5','1.0','1.5','2'})
% clim([0 2])
% caxis=([0 2])
% title("102 정선 염분 표준편차 ( psu) 1980 - 2019",'FontSize',20,'FontWeight','bold')
% print('std102S','-depsc','-tiff')
% 
% colormap("parula")
% colorbar
% clim([0 1])
% caxis([0 1])
% title("102 정선 DO 표준편차 ( mL/L) 1980 - 2019",'FontSize',20,'FontWeight','bold')
% print('std102DO','-depsc','-tiff')
%% plot vertical section maximum
% max102_04 = max(vq,[],2,'omitnan');
% figure(3)
% plot(max102_04,xq,'ro')
% set(gca,'ydir','rev')

%% plot vertical section minimum
% min102_04 = min(vq,[],2,'omitnan');
% figure(4)
% plot(min102_04,xq,'bo')
% set(gca,'ydir','rev')

%% plot linear trend (trend of what? section? horizontal surface? mass?)


%% plot number of missing values


%% plot number of QC==2


%% plot number of QC==3





