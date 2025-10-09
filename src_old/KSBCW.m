%% KSBCW Long term variability
%% Temperature
% clear data

clc; close all; clear
% load variables

% load mat file
load temperature.mat; % T
lev = [0 10 20 30 50 75 100 125 150 200 250 300 400 500]; % Standard depth level
depth_max_208 = [70 130 150 143 140 135 130 130 130]; % 01, 02, 03, 04, 05, 06, 07, 08, 09
depth_max_207 = [84 103 115 227 160]; % 01, 02, 03, 04, 05
% calculate distance between the stations

lats = [35.475 35.405 35.2967 35.185 35.08 34.97 34.8667 34.755 34.6483]; % 북위
lons = [129.455 129.5633 129.715 129.8783 130.033 130.1983 130.3483 130.51 130.675]; % 동경
d = nan(1,9);
d(1) = 0;
for i=1:8
    d_temp = sqrt(lons(i+1)^2 + lats(i+1)^2) - sqrt(lons(i)^2 +lats(i)^2);
    d(i+1) = d_temp;
    clearvars d_temp
end
% Line 208 station 01 설정

target_line = 208;
station_num = unique(T(T(:,1)==target_line,2));
t_5 = [];
for l=1:length(station_num)
     target_station = station_num(l);
    i_line = T(:,1) == target_line;
    i_station = T(:,2) == target_station;
    depth_max = depth_max_208(target_station);
    depths = T(i_line&i_station,3);
    dates = T(i_line&i_station,4:6);
    c_dates = categorical(datetime(dates));
    t_ob = T(i_line&i_station,7);
    qc_flags = T(i_line&i_station,8);
    clearvars i_line i_station
    if l==6
        depths = [depths(1:842);depths(851:end)];
        dates = [dates(1:842,:);dates(851:end,:)];
        c_dates = categorical(datetime(dates));
        t_ob = [t_ob(1:842);t_ob(851:end)];
        qc_flags = [qc_flags(1:842);qc_flags(851:end)];
    elseif l==8
        depths = [depths(1:200);depths(209:end)];
        dates = [dates(1:200,:);dates(209:end,:)];
        c_dates = categorical(datetime(dates));
        t_ob = [t_ob(1:200);t_ob(209:end)];
        qc_flags = [qc_flags(1:200);qc_flags(209:end)];
    elseif l==9
        depths = [depths(1:200);depths(209:end)];
        dates = [dates(1:200,:);dates(209:end,:)];
        c_dates = categorical(datetime(dates));
        t_ob = [t_ob(1:200);t_ob(209:end)];
        qc_flags = [qc_flags(1:200);qc_flags(209:end)];
    end
% 1 m 간격 깊이 반영하기

    t_1 = nan(depth_max+1,480);
    qc_reshape = nan(size(t_1));
    for i=0:500
        i_depth = depths==i;
        for j=1:40
            i_year = dates(:,1) == 1979+j;
            for k=1:12
                i_month= dates(:,2) ==k;
                if qc_flags(i_depth&i_year&i_month)==4;
                    t_ob(i_depth&i_year&i_month)=NaN;
                end
                if ~isnan(t_ob(i_depth&i_year&i_month))
                    t_1(i+1,k+12*(j-1)) = unique(t_ob(i_depth&i_year&i_month)); % observation data depth*time
                    qc_reshape(i+1,k+12*(j-1)) = unique(qc_flags(i_depth&i_year&i_month));
                end
            end
        end
    end
    clearvars i j k i_depth i_year i_month
    t_1(t_1<0.02)=NaN;
% QC Flag counts

    [num_qc1,edges] = histcounts(c_dates(qc_flags==1),categories(c_dates));
    num_qc2 = histcounts(c_dates(qc_flags==2),categories(c_dates));
    percent_qc = num_qc2./(num_qc1+num_qc2)*100;
    qc_others_t(l) = {c_dates(qc_flags==3|qc_flags==4|qc_flags==5)};
    [~,col] = find(qc_reshape==4);
    figure;
    bar(percent_qc(percent_qc~=0));
    xticks(1:length(percent_qc(percent_qc~=0)))
    xticklabels(edges(percent_qc~=0));
    title(strcat("Flag percent of line 208 station",string(station_num(l))),"FontSize",20,"FontWeight","bold")
    ylabel("Percent ( %)","FontSize",18,"FontWeight","bold")
    box("off")
    grid on
    exportgraphics(gca,strcat("flag_count",string(l),".tiff"),"ContentType","auto","BackgroundColor","white")
    close(gcf)
    t_2 = reshape(t_1(:,:),[],12,40); % observation data depth*month*year
% Depth interpolation

    t_3 = nan(length(lev),480);
    for i=1:480
        if ~isnan(t_1(1,i))
            for j=1:nnz(lev<depth_max)
                F = griddedInterpolant(find(~isnan(t_1(1:depth_max+1,i))),t_1(~isnan(t_1(1:depth_max+1,i)),i),"linear");
                vq = F(lev(lev<depth_max)+1);
                t_3(j,i) = vq(j);
                clearvars F vq
            end
        end
    end
    clearvars i j
    t_3(t_3<0|t_3>33)=NaN;
% Depth interpolation check
    figure;
        for i=1:480
            if nnz(col==i)>=1
                plot(t_1(~isnan(t_1(1:depth_max+1,i)),i),find(~isnan(t_1(1:depth_max+1,i)))-1,"LineStyle","none","LineWidth",1.5,"Marker","o","MarkerSize",10,"MarkerEdgeColor","r")
                hold on
                plot(t_3(:,i),lev,"LineStyle","none","LineWidth",1.5,"Marker","x","MarkerSize",10,"MarkerEdgeColor","b")
                plot(t_3(:,i),lev,"-k","LineWidth",1.0)
                hold off
                xtickformat("%2.1f");
                legend("Raw","Interpolation","Location","southeast")
                xlabel("Temperature ( \circC)","FontSize",14,"FontWeight","bold")
                ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
                set(gca,"YDir","reverse")
                grid on
                if mod(i,12)~=0
                    title(strcat("Station",num2str(l,'%02.0f'),"\_",num2str(1980+floor(i/12)),"-",num2str(mod(i,12),"%02.0f")),"FontSize",20,"FontWeight","bold");
                    exportgraphics(gca,strcat("Depth_interpolation_","Station",num2str(l,'%02.0f'),"_",num2str(1980+floor(i/12),"%4.0f"),"-",num2str(mod(i,12),"%02.0f"),"_temperature",".tiff"),"ContentType","auto","BackgroundColor","white")
                    close(gcf)
                else
                    title(strcat("Station",num2str(l,'%02.0f'),"\_",num2str(1979+floor(i/12)),"-",num2str(12,"%02.0f")),"FontSize",20,"FontWeight","bold");
                    exportgraphics(gca,strcat("Depth_interpolation_","Station",num2str(l,'%02.0f'),"_",num2str(1979+floor(i/12),"%4.0f"),"-",num2str(12,"%2.0f"),"_temperature",".tiff"),"ContentType","auto","BackgroundColor","white")
                    close(gcf)
                end
            end
        end

% Time interpolation

    [row,col]=size(t_3);
    t_4=nan(row,col/2);
    for i=1:row
        for j=1:find(~isnan(t_3(i,:)),1,'last')
            F=griddedInterpolant(find(~isnan(t_3(i,:))),t_3(i,~isnan(t_3(i,:))),"linear");
            vq = F(1:find(~isnan(t_3(i,:)),1,'last'));
            if rem(j,2)==0
               t_4(i,j/2)=vq(j);
            end
            clearvars F vq
        end
    end
    clearvars row col i j
% time interpolation check

    for i=1:length(lev)
        if ~isnan(t_4(i,2))
            figure;
            plot(t_1(lev(i)+1,2:2:480),"LineStyle","none","LineWidth",1.5,"Marker","o","MarkerSize",10,"MarkerEdgeColor","r")
            hold on
            plot(1:240,t_4(i,:),"LineStyle","none","LineWidth",1.5,"Marker","x","MarkerSize",10,"MarkerEdgeColor","b")
            plot(1:240,t_4(i,:),'-k',"LineWidth",1.0)
            hold off
            xlim([0 240])
            ylim([-5 40])
            xticks(linspace(0,240,5))
            xticklabels(["1980","1990","2000","2010","2020"])
            xlabel("Year","FontSize",14,"FontWeight","bold")
            ylabel("Temperature (\circC)","FontSize",14,"FontWeight","bold");
            legend('Raw','Interpolation',"Location","northeast")
            title(strcat(string(lev(i)),"m"),"FontSize",20,"FontWeight","bold");
            grid on
            exportgraphics(gca,strcat("time_interpolation_Temperature",string(l),string(lev(i)),".tiff"),"ContentType","auto","BackgroundColor","white")
            close(gcf)
        end
    end
    clearvars i
    t_5 = [t_5,t_4];
% Climatology

    clima = squeeze(mean(reshape(t_4,14,6,[]),3,"omitnan"));
    figure;
    plot(clima',"LineWidth",1.0)
    for i=1:nnz(lev<depth_max)
        legendlabels(i) = strcat(string(lev(i))," m");
    end
    clearvars i
    legend(legendlabels)
    grid on
    ylim([0 30])
    ylabel("Temperature ( \circC)","FontSize",14,"FontWeight","bold")
    xlim([1 6])
    xticks(1:1:6)
    xticklabels(["Feburary","April","June","August","October","December"])
    xtickangle(45)
    title("Climatology","FontSize",20,"FontWeight","bold")
    exportgraphics(gca,strcat("climatology",string(l),".tiff"),"ContentType","auto","BackgroundColor","white")
    close(gcf)
end

% vertical mean section

t_6 = reshape(t_5,14,6,40,9);
t_7 = squeeze(mean(t_6(:,:,1:19,:),3,"omitnan")); % 1980-1998, all station
t_8 = squeeze(mean(t_7,2,"omitnan"));
t_9 = squeeze(mean(t_6(:,:,:,1:4),3,"omitnan")); % 1-4 station, all period
t_10 = squeeze(mean(t_9,2,"omitnan"));
t_11 = squeeze(mean(t_6(:,:,20:end,1:4),3,"omitnan")); % 1-4 station, 1999-2019
t_12 = squeeze(mean(t_11,2,"omitnan"));

tmin = 0; tdif = 1; tmax = 30;

figure; % 1980-1998 mean temperature
contourf(cumsum(d),lev,t_8,tmin:tdif/3:tmax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d),lev,t_8,tmin:3*tdif:tmax,"LineStyle","-","LineColor","k","ShowText","on")
contour(cumsum(d),lev,t_8,[10 10],"LineStyle","-","LineColor","k","LineWidth",2.0,"ShowText","on")
grid on
hold off
xticks(cumsum(d))
xticklabels(["01","02","03","04","05","06","07","08","09"])
colorbar("Ticks",tmin:3*tdif:tmax);
caxis([0 21])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean Temperature (\circC) 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_mean_temperature.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-2019 mean temperature
contourf(cumsum(d(1:4)),lev,t_10,tmin:tdif/3:tmax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,t_10,tmin:3*tdif:tmax,"LineStyle","-","LineColor","k","ShowText","on")
contour(cumsum(d(1:4)),lev,t_10,[10 10],"LineStyle","-","LineColor","k","LineWidth",2.0,"ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",tmin:3*tdif:tmax);
caxis([0 21])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean Temperature (\circC) 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_mean_temperature.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 mean temperature
contourf(cumsum(d(1:4)),lev,t_12,tmin:tdif/3:tmax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,t_12,tmin:3*tdif:tmax,"LineStyle","-","LineColor","k","ShowText","on")
contour(cumsum(d(1:4)),lev,t_12,[10 10],"LineStyle","-","LineColor","k","LineWidth",2.0,"ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",tmin:3*tdif:tmax);
caxis([0 21])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean Temperature (\circC) 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_mean_temperature.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

% vertical std section

std_8 = squeeze(std(t_7,0,2,"omitnan"));
std_10 = squeeze(std(t_9,0,2,"omitnan"));
std_12 = squeeze(std(t_11,0,2,"omitnan"));

figure; % 1980-2019 temperature standard deviation
contourf(cumsum(d(1:4)),lev,std_10,0:.05:5,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,std_10,0:1:5,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:1:5);
caxis([0 5])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Temperature Std (\circC) 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_std_temperature.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 temperature standard deviation
contourf(cumsum(d),lev,std_8,0:.05:5,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d),lev,std_8,0:1:5,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d))
xticklabels(["01","02","03","04","05","06","07","08","09"])
colorbar("Ticks",0:1:5);
caxis([0 5])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Temperature Std (\circC) 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_std_temperature.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 temperature standard deviation
contourf(cumsum(d(1:4)),lev,std_12,0:.05:5,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,std_12,0:1:5,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:1:5);
caxis([0 5])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Temperature Std (\circC) 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_std_temperature.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)
%% Salinity

load Salt; % S
s_5 = [];
for l=1:length(station_num)
     target_station = station_num(l);
    i_line = S(:,1) == target_line;
    i_station = S(:,2) == target_station;
    depth_max = depth_max_208(target_station);
    depths = S(i_line&i_station,3);
    dates = S(i_line&i_station,4:6);
    c_dates = categorical(datetime(dates));
    s_ob = S(i_line&i_station,7);
    qc_flags = S(i_line&i_station,8);
    clearvars i_line i_station
     if l==6
         depths = [depths(1:842);depths(851:end)];
         dates = [dates(1:842,:);dates(851:end,:)];
         c_dates = categorical(datetime(dates));
         s_ob = [s_ob(1:842);s_ob(851:end)];
         qc_flags = [qc_flags(1:842);qc_flags(851:end)];
     elseif l==8
         depths = [depths(1:200);depths(209:end)];
         dates = [dates(1:200,:);dates(209:end,:)];
         c_dates = categorical(datetime(dates));
         s_ob = [s_ob(1:200);s_ob(209:end)];
         qc_flags = [qc_flags(1:200);qc_flags(209:end)];
     elseif l==9
         depths = [depths(1:200);depths(209:end)];
         dates = [dates(1:200,:);dates(209:end,:)];
         c_dates = categorical(datetime(dates));
         s_ob = [s_ob(1:200);s_ob(209:end)];
         qc_flags = [qc_flags(1:200);qc_flags(209:end)];
     end
%      s_ob(s_ob>40) = NaN;
% 1 m 간격 깊이 반영하기

    s_1 = nan(depth_max+1,480);
    qc_reshape = nan(size(s_1));
    for i=0:500
        i_depth = depths==i;
        for j=1:40
            i_year = dates(:,1) == 1979+j;
            for k=1:12
                i_month= dates(:,2) ==k;
                if ~isnan(s_ob(i_depth&i_year&i_month))
                    if qc_flags(i_depth&i_year&i_month)==4
                        s_ob(i_depth&i_year&i_month)=NaN;
                    end
                    s_1(i+1,k+12*(j-1)) = unique(s_ob(i_depth&i_year&i_month)); % observation data depth*time
                    qc_reshape(i+1,k+12*(j-1)) = unique(qc_flags(i_depth&i_year&i_month));
                end
            end
        end
    end
    clearvars i j k i_depth i_year i_month
    s_1(s_1<28)=NaN;
% QC Flag counts

    [num_qc1,edges] = histcounts(c_dates(qc_flags==1),categories(c_dates));
    num_qc2 = histcounts(c_dates(qc_flags==2),categories(c_dates));
    percent_qc = num_qc2./(num_qc1+num_qc2)*100;
    qc_others_s(l) = {c_dates(qc_flags==3|qc_flags==4|qc_flags==5)};
    [~,col] = find(qc_reshape==4);
    figure; % draw flag
    bar(percent_qc(percent_qc~=0));
    xticks(1:length(percent_qc(percent_qc~=0)))
    xticklabels(edges(percent_qc~=0));
    title(strcat("Flag percent of line 208 station",string(station_num(l))),"FontSize",20,"FontWeight","bold")
    ylabel("Percent ( %)","FontSize",18,"FontWeight","bold")
    box("off")
    grid on
    exportgraphics(gca,strcat("salt_flag_count",string(l),".tiff"),"ContentType","auto","BackgroundColor","white")
    close(gcf)
    s_2 = reshape(s_1(:,:),[],12,40); % observation data depth*month*year

% Depth interpolation

    s_3 = nan(length(lev),480);
    for i=1:480
        if ~isnan(s_1(1,i))
            for j=1:nnz(lev<depth_max)
                F = griddedInterpolant(find(~isnan(s_1(1:depth_max+1,i))),s_1(~isnan(s_1(1:depth_max+1,i)),i),"linear");
                vq = F(lev(lev<depth_max)+1);
                s_3(j,i) = vq(j);
                clearvars F vq
            end
        elseif isnan(s_1(1,i))&&~isnan(s_1(11,i))
            for j=1:nnz(lev<depth_max)
                F = griddedInterpolant(find(~isnan(s_1(1:depth_max+1,i))),s_1(~isnan(s_1(1:depth_max+1,i)),i),"linear");
                vq = F(lev(lev<depth_max)+1);
                s_3(j,i) = vq(j);
                clearvars F vq
            end
        end
    end
    clearvars i j
    s_3(s_3>35)=nan;

%Depth interpolation check
    figure; % depth interpolation of salinity
    for i=1:480
        if nnz(col==i)>=1
            plot(s_1(~isnan(s_1(1:depth_max+1,i)),i),find(~isnan(s_1(1:depth_max+1,i)))-1,"LineStyle","none","LineWidth",1.5,"Marker","o","MarkerSize",10,"MarkerEdgeColor","r")
            hold on
            plot(s_3(:,i),lev,"LineStyle","none","LineWidth",1.5,"Marker","x","MarkerSize",10,"MarkerEdgeColor","b")
            plot(s_3(:,i),lev,"-k","LineWidth",1.0)
            hold off
            xtickformat("%2.3f");
            legend("Raw","Interpolation","Location","southeast")
            xlabel("Salinity (psu)","FontSize",14,"FontWeight","bold")
            ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
            set(gca,"YDir","reverse")
            grid on
            if mod(i,12)~=0
                title(strcat("Station",num2str(l,'%02.0f'),"\_",num2str(1980+floor(i/12)),"-",num2str(mod(i,12),"%02.0f")),"FontSize",20,"FontWeight","bold");
                exportgraphics(gca,strcat("Depth_interpolation_","Station",num2str(l,'%02.0f'),"_",num2str(1980+floor(i/12),"%4.0f"),"-",num2str(mod(i,12),"%02.0f"),"_Salinity",".tiff"),"ContentType","auto","BackgroundColor","white")
                close(gcf)
            else
                title(strcat("Station",num2str(l,'%02.0f'),"\_",num2str(1979+floor(i/12)),"-",num2str(12,"%02.0f")),"FontSize",20,"FontWeight","bold");
                exportgraphics(gca,strcat("Depth_interpolation_","Station",num2str(l,'%02.0f'),"_",num2str(1979+floor(i/12),"%4.0f"),"-",num2str(12,"%2.0f"),"_Salinity",".tiff"),"ContentType","auto","BackgroundColor","white")
                close(gcf)
            end
        end
    end
% Time interpolation

    [row,col]=size(s_3);
    s_4=nan(row,col/2);
    for i=1:row
        for j=1:find(~isnan(s_3(i,:)),1,'last')
            F=griddedInterpolant(find(~isnan(s_3(i,:))),s_3(i,~isnan(s_3(i,:))),"linear");
            vq = F(1:find(~isnan(s_3(i,:)),1,'last'));
            if rem(j,2)==0
               s_4(i,j/2)=vq(j);
            end
            clearvars F vq
        end
    end
    clearvars row col i j
% time interpolation check
    
    for i=1:length(lev)
        if ~isnan(s_4(i,2))
            figure; % time interpolation of salinity
            plot(s_1(lev(i)+1,2:2:480),"LineStyle","none","LineWidth",1.5,"Marker","o","MarkerSize",10,"MarkerEdgeColor","r")
            hold on
            plot(1:240,s_4(i,:),"LineStyle","none","LineWidth",1.5,"Marker","x","MarkerSize",10,"MarkerEdgeColor","b")
            plot(1:240,s_4(i,:),'-k',"LineWidth",1.0)
            hold off
            xlim([0 240])
            ylim([28 36])
            xticks(linspace(0,240,5))
            xticklabels(["1980","1990","2000","2010","2020"])
            xlabel("Year","FontSize",14,"FontWeight","bold")
            ylabel("Salinity","FontSize",14,"FontWeight","bold");
            legend('Raw','Interpolation',"Location","northeast")
            title(strcat(string(lev(i)),"m"),"FontSize",20,"FontWeight","bold");
            grid on
            exportgraphics(gca,strcat("time_interpolation_salt",string(l),string(lev(i)),".tiff"),"ContentType","auto","BackgroundColor","white")
            close(gcf)
        end
    end
    clearvars i
    s_5 = [s_5,s_4];
% Climatology

    clima_s = squeeze(mean(reshape(s_4,14,6,[]),3,"omitnan"));
    figure; % salinity climatology
    plot(clima_s',"LineWidth",1.0)
    for i=1:nnz(lev<depth_max)
        legendlabels(i) = strcat(string(lev(i))," m");
    end
    clearvars i
    legend(legendlabels)
    grid on
    ylim([28 36])
    ylabel("Salinity","FontSize",14,"FontWeight","bold")
    xlim([1 6])
    xticks(1:1:6)
    xticklabels(["Feburary","April","June","August","October","December"])
    xtickangle(45)
    title("Climatology S","FontSize",20,"FontWeight","bold")
    exportgraphics(gca,strcat("climatology_s",string(l),".tiff"),"ContentType","auto","BackgroundColor","white")
    close(gcf)
end

    qc_s_mat = [];
    for i=1:9
        qc_s_mat = vertcat(qc_s_mat,qc_others_s{1,i});
    end

% vertical mean section

s_6 = reshape(s_5,14,6,40,9);
s_7 = squeeze(mean(s_6(:,:,1:19,:),3,"omitnan"));
s_8 = squeeze(mean(s_7,2,"omitnan"));
s_9 = squeeze(mean(s_6(:,:,:,1:4),3,"omitnan"));
s_10 = squeeze(mean(s_9,2,"omitnan"));
s_11 = squeeze(mean(s_6(:,:,20:end,1:4),3,"omitnan"));
s_12 = squeeze(mean(s_11,2,"omitnan"));
smin = 33; sdif=1;smax=35;

figure; % 1980-1998 mean salinity
contourf(cumsum(d),lev,s_8,smin:sdif/100:smax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d),lev,s_8,smin:sdif/10:smax,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d))
xticklabels(["01","02","03","04","05","06","07","08","09"])
colorbar("Ticks",smin:sdif/10:smax);
 caxis([33.6 34.6])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean Salinity 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_mean_salinity.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-2018 mean salinity
contourf(cumsum(d(1:4)),lev,s_10,smin:sdif/100:smax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,s_10,smin:sdif/10:smax,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",smin:sdif/10:smax);
caxis([33.6 34.6])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean Salinity 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_mean_salinity.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2018 mean salinity
contourf(cumsum(d(1:4)),lev,s_12,smin:sdif/100:smax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,s_12,smin:sdif/10:smax,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",smin:sdif/10:smax);
caxis([33.6 34.6])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean Salinity 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_mean_salinity.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

% vertical std section

sstd_7 = squeeze(std(s_6(:,:,1:19,:),0,3,"omitnan"));
sstd_8 = squeeze(std(s_7,0,2,"omitnan"));
sstd_9 = squeeze(std(s_6(:,:,:,1:4),0,3,"omitnan"));
sstd_10 = squeeze(std(s_9,0,2,"omitnan"));
sstd_11 = squeeze(std(s_6(:,:,20:end,1:4),0,3,"omitnan"));
sstd_12 = squeeze(std(s_11,0,2,"omitnan"));


figure; % 1980-1998 salinity standard deviation
contourf(cumsum(d),lev,sstd_8,0:.01:2,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d),lev,sstd_8,0:.2:2,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d))
xticklabels(["01","02","03","04","05","06","07","08","09"])
colorbar("Ticks",0:.2:2)
caxis([0 1])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Salinity std 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_std_salinity.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-2019 salinity standard deviation
contourf(cumsum(d(1:4)),lev,sstd_10,0:.01:2,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,sstd_10,0:.2:2,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:.2:2)
caxis([0 1])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Salinity std 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_std_salinity.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 salinity standard deviation
contourf(cumsum(d(1:4)),lev,sstd_12,0:.01:2,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,sstd_12,0:.2:2,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:.2:2)
caxis([0 1])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Salinity std 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_std_salinity.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)
%% CT-SA

% pressure
p = gsw_p_from_z(repmat(-lev',1,length(lats)),lats);

% absolute salinity
sa = nan(14,240,9);
s_11 = reshape(s_5,14,240,9);
for i=1:240
    sa(:,i,:) = gsw_SA_from_SP(squeeze(s_11(:,i,:)),p,lons,lats);
end

% conservative temperature
ct = nan(14,240,9);
t_11 = reshape(t_5,14,240,9);
for i=1:240
    ct(:,i,:) = gsw_CT_from_t(squeeze(sa(:,i,:)),squeeze(t_11(:,i,:)),p);
end

% draw

% Mean
figure; % 1980-2019 conservative temperature mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(ct(:,:,1:4),2,"omitnan")),0:.1:21,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(mean(ct(:,:,1:4),2,"omitnan")),0:3:21,"LineStyle","-","LineColor","k","ShowText","on")
contour(cumsum(d(1:4)),lev,squeeze(mean(ct(:,:,1:4),2,"omitnan")),[10 10],"LineStyle","-","LineColor","k","LineWidth",2.0,"ShowText","on")
set(gca,"ydir","reverse")
colormap("jet");
ylim([0 125])
clim([0 21])
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("ticks",0:3:21)
caxis([0 21])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean CT (\circC) 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_mean_CT.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 conservative temperature mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(ct(:,1:114,1:4),2,"omitnan")),0:.1:21,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(mean(ct(:,1:114,1:4),2,"omitnan")),0:3:21,"LineStyle","-","LineColor","k","ShowText","on")
contour(cumsum(d(1:4)),lev,squeeze(mean(ct(:,1:114,1:4),2,"omitnan")),[10 10],"LineStyle","-","LineColor","k","LineWidth",2.0,"ShowText","on")
set(gca,"ydir","reverse")
colormap("jet");
ylim([0 125])
clim([0 21])
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("ticks",0:3:21)
caxis([0 21])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean CT (\circC) 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_mean_CT.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 conservative temperature mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(ct(:,115:end,1:4),2,"omitnan")),0:.1:21,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(mean(ct(:,115:end,1:4),2,"omitnan")),0:3:21,"LineStyle","-","LineColor","k","ShowText","on")
contour(cumsum(d(1:4)),lev,squeeze(mean(ct(:,115:end,1:4),2,"omitnan")),[10 10],"LineStyle","-","LineColor","k","LineWidth",2.0,"ShowText","on")
set(gca,"ydir","reverse")
colormap("jet");
ylim([0 125])
clim([0 21])
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("ticks",0:3:21)
caxis([0 21])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean CT (\circC) 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_mean_CT.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-2019 absolute salinity mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(sa(:,:,1:4),2,"omitnan")),smin:.01:smax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,squeeze(mean(sa(:,:,1:4),2,"omitnan")),smin:sdif/10:smax,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",smin:sdif/10:smax);
caxis([33.6 34.6])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean SA 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_mean_SA.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 absolute salinity mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(sa(:,1:114,1:4),2,"omitnan")),smin:.01:smax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,squeeze(mean(sa(:,1:114,1:4),2,"omitnan")),smin:sdif/10:smax,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",smin:sdif/10:smax);
caxis([33.6 34.6])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean SA 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_mean_SA.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999 - 2019 absolute salinity mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(sa(:,114:end,1:4),2,"omitnan")),smin:.01:smax,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,squeeze(mean(sa(:,114:end,1:4),2,"omitnan")),smin:sdif/10:smax,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",smin:sdif/10:smax);
caxis([33.6 34.6])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Mean SA 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_mean_SA.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

% density
rho = nan(14,240,9);
for i=1:240
    rho(:,i,:) = gsw_rho(squeeze(sa(:,i,:)),squeeze(ct(:,i,:)),p);
end

figure; % 1980-2019 density mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(rho(:,:,1:4),2,"omitnan")-1000),24:.1:28,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(mean(rho(:,:,1:4),2,"omitnan")-1000),24:1:28,"LineStyle","-","LineColor","k","ShowText","on")
hold off
colormap("jet")
grid on
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar('Ticks',24:1:28)
caxis([24 28])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("\sigma (kg/m^3) 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_rho.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 density mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(rho(:,1:114,1:4),2,"omitnan")-1000),24:.1:28,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(mean(rho(:,1:114,1:4),2,"omitnan")-1000),24:1:28,"LineStyle","-","LineColor","k","ShowText","on")
hold off
colormap("jet")
grid on
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar('Ticks',24:1:28)
caxis([24 28])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("\sigma (kg/m^3) 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_rho.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 density mean
contourf(cumsum(d(1:4)),lev,squeeze(mean(rho(:,115:end,1:4),2,"omitnan")-1000),23:.1:28,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(mean(rho(:,115:end,1:4),2,"omitnan")-1000),23:1:28,"LineStyle","-","LineColor","k","ShowText","on")
hold off
colormap("jet")
grid on
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar('Ticks',24:1:28)
caxis([24 28])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("\sigma (kg/m^3) 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_rho.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

% Standard deviation
figure; % 1980-2019 conservative temperature standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(ct(:,:,1:4),0,2,"omitnan")),0:.05:5,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(std(ct(:,:,1:4),0,2,"omitnan")),0:1:5,"LineStyle","-","LineColor","k","ShowText","on")
set(gca,"ydir","reverse")
colormap("jet");
ylim([0 125])
clim([0 5])
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("ticks",0:1:5)
caxis([0 5])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Std CT (\circC) 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_std_CT.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 conservative temperature standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(ct(:,1:114,1:4),0,2,"omitnan")),0:.05:5,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(std(ct(:,1:114,1:4),0,2,"omitnan")),0:1:5,"LineStyle","-","LineColor","k","ShowText","on")
set(gca,"ydir","reverse")
colormap("jet");
ylim([0 125])
clim([0 5])
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("ticks",0:1:5)
caxis([0 5])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Std CT (\circC) 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_std_CT.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 conservative temperature mean
contourf(cumsum(d(1:4)),lev,squeeze(std(ct(:,115:end,1:4),0,2,"omitnan")),0:.05:5,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(std(ct(:,115:end,1:4),0,2,"omitnan")),0:1:5,"LineStyle","-","LineColor","k","ShowText","on")
set(gca,"ydir","reverse")
colormap("jet");
ylim([0 125])
clim([0 5])
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("ticks",0:1:5)
caxis([0 5])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Std CT (\circC) 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_std_CT.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-2019 absolute salinity standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(sa(:,:,1:4),0,2,"omitnan")),0:.01:2,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,squeeze(std(sa(:,:,1:4),0,2,"omitnan")),0:.2:2,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:.2:1);
caxis([0 1])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Std SA 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_std_SA.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 absolute salinity standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(sa(:,1:114,1:4),0,2,"omitnan")),0:.01:2,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,squeeze(std(sa(:,1:114,1:4),0,2,"omitnan")),0:.2:2,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:.2:1);
caxis([0 1])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Std SA 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_std_SA.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999 - 2019 absolute salinity standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(sa(:,114:end,1:4),0,2,"omitnan")),0:.01:2,"LineStyle","none")
hold on
colormap("jet")
contour(cumsum(d(1:4)),lev,squeeze(std(sa(:,114:end,1:4),0,2,"omitnan")),0:.2:2,"LineStyle","-","LineColor","k","ShowText","on")
grid on
hold off
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar("Ticks",0:.2:1);
caxis([0 1])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("Std SA 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_std_SA.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

% density
rho = nan(14,240,9);
for i=1:240
    rho(:,i,:) = gsw_rho(squeeze(sa(:,i,:)),squeeze(ct(:,i,:)),p);
end

figure; % 1980-2019 density standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(rho(:,:,1:4),0,2,"omitnan")),0:.02:2,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(std(rho(:,:,1:4),0,2,"omitnan")),0:.4:2,"LineStyle","-","LineColor","k","ShowText","on")
hold off
colormap("jet")
grid on
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar('Ticks',0:.4:2)
caxis([0 2])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("\sigma (kg/m^3) std 1980 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-2019_vertical_rho_std.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1980-1998 density standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(rho(:,1:114,1:4),0,2,"omitnan")),0:.02:2,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(std(rho(:,1:114,1:4),0,2,"omitnan")),0:.4:2,"LineStyle","-","LineColor","k","ShowText","on")
hold off
colormap("jet")
grid on
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar('Ticks',0:.4:2)
caxis([0 2])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("\sigma (kg/m^3) std 1980 - 1998","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1980-1998_vertical_rho_std.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

figure; % 1999-2019 density standard deviation
contourf(cumsum(d(1:4)),lev,squeeze(std(rho(:,115:end,1:4),0,2,"omitnan")),0:.02:2,"LineStyle","none")
hold on
contour(cumsum(d(1:4)),lev,squeeze(std(rho(:,115:end,1:4),0,2,"omitnan")),0:.4:2,"LineStyle","-","LineColor","k","ShowText","on")
hold off
colormap("jet")
grid on
xticks(cumsum(d(1:4)))
xticklabels(["01","02","03","04"])
colorbar('Ticks',0:.4:2)
caxis([0 2])
ylim([0 125])
yticks(lev)
ylabel("Depth (m)","FontSize",14,"FontWeight","bold")
set(gca,"ydir","reverse")
xlabel("Station","FontSize",14,"FontWeight","bold")
title("\sigma (kg/m^3) std 1999 - 2019","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"1999-2019_vertical_rho_std.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)
%% T-S diagram

[Tmesh, Smesh] = meshgrid(tmin:(tmax-tmin)/13:tmax,smin:(smax-smin)/13:smax);
[CTmesh, SAmesh] = meshgrid(squeeze(mean(ct(:,:,1:4),2,"omitnan")),squeeze(mean(sa(:,:,1:4),2,"omitnan")));
rho_mesh = gsw_rho(Smesh,Tmesh,p(:,1));
figure;
contour(Smesh,Tmesh,rho_mesh,"LineColor","k","ShowText","on")
hold on
plot(squeeze(mean(sa(:,:,1:4),2,"omitnan")),squeeze(mean(ct(:,:,1:4),2,"omitnan")),".k","MarkerSize",15);
plot(squeeze(mean(sa(:,1:114,1:4),2,"omitnan")),squeeze(mean(ct(:,1:114,1:4),2,"omitnan")),".b","MarkerSize",15);
plot(squeeze(mean(sa(:,115:end,1:4),2,"omitnan")),squeeze(mean(ct(:,115:end,1:4),2,"omitnan")),".r","MarkerSize",15);
yline(10,"k-","LineWidth",1.5)
xline(34.5,"k-","LineWidth",1.5)
hold off
xtickformat("%3.1f");
xlabel("Absolute Salinity","FontSize",14,"FontWeight","bold")
ylabel("Conservative Temperature (°C)","FontSize",14,"FontWeight","bold")
title("T-S Diagram","FontSize",20,"FontWeight","bold")
exportgraphics(gca,"TSdiagram.tiff","ContentType","auto","BackgroundColor","white")
close(gcf)

%% EOF
