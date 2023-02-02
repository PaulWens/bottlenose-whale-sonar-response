%% Plot Acoustic detections in JM1 and JM5 recordings and ECDF of JM5 pre-exposure
clear all; close all;

dt = [5 8]; % detector thresholds
col = [.6 .6 .6; 0 0 0]; % grey and black
yoffset = 0.2;
tson16 = datenum({'18-Jun-2016 12:16:00','18-Jun-2016 12:51:25'});
tson151 = datenum({'20-Jun-2015 15:13:00','20-Jun-15 15:28:00'}); 
tson152 = datenum({'29-Jun-2015 02:48:00','29-Jun-15 03:03:00'}); 

%% DETECTIONS TIME SERIES 2015
figure; subplot(2,2,1)
[tmp,~] = xlsread('NBWdetections.xlsx','JM1');
time = tmp(:,1); snr = tmp(:,2);
for i=1:2
nbwDet = snr>dt(i);
ts = floor(time(1));
    for j=1:20
        ind = time>=floor(ts)+(j-1) & time<floor(ts)+j & nbwDet;
        timeofday = time(ind)-floor(min(time(ind)));
        tmp = datevec(min(time(ind)));
        if ~isempty(timeofday)
            hold on; plot(timeofday,nbwDet(ind)*tmp(3)-yoffset*(i-1),'s','Color',col(i,:),'MarkerFaceColor',col(i,:),'MarkerSize',4); hold off
        end
    end
end
datetick('x',15)
ylabel('Day of month'); xlabel('Time of day (UTC)'); 
tmp1=datevec(time(1)); tmp2=datevec(time(end)); 
set(gca,'YLim',[tmp1(3)-1 tmp2(3)+1],'YTick',16:2:30)

% Plot start and end of recording
hold on
tmp1=datevec(time([1 end])); tmp2=(time([1 end])-floor(time([1 end])));
plot(tmp2,tmp1(:,3),'x','LineWidth',2,'MarkerSize',8,'Color',col(1,:));

% Plot sonar periods
tmp1=datevec(tson151); tmp2=(tson151-floor(tson151));
plot(tmp2,tmp1(:,3)+yoffset,'r-','LineWidth',6); 
tmp1=datevec(tson152); tmp2=(tson152-floor(tson152));
plot(tmp2,tmp1(:,3)+yoffset,'r-','LineWidth',6); 
hold off

%% DETECTIONS TIME SERIES 2016
subplot(2,2,2)
offset = [-0.1 0 0.1];
[tmp,~] = xlsread('NBWdetections.xlsx','JM5');
time = tmp(:,1); snr = tmp(:,2);
for i=1:2
    nbwDet = snr>dt(i);
    ts = floor(time(1));
    for j=1:20
        ind = time>=floor(ts)+(j-1) & time<floor(ts)+j & nbwDet;
        timeofday = time(ind)-floor(min(time(ind)));
        tmp = datevec(min(time(ind)));
        if ~isempty(timeofday)
            hold on; plot(timeofday,nbwDet(ind)*tmp(3)-yoffset*(i-1),'s','Color',col(i,:),'MarkerFaceColor',col(i,:),'MarkerSize',4); hold off
        end
    end
end
datetick('x',15)
ylabel('Day of month'); xlabel('Time of day (UTC)');
tmp1=datevec(time(1)); tmp2=datevec(time(end)); 
set(gca,'YLim',[tmp1(3)-1 tmp2(3)+1])

% Plot start and end of recording
hold on
tmp1=datevec(time([1 end])); tmp2=(time([1 end])-floor(time([1 end])));
plot(tmp2,tmp1(:,3),'x','LineWidth',2,'MarkerSize',8,'Color',col(1,:));

% Plot sonar period
tmp1=datevec(tson16); tmp2=(tson16-floor(tson16));
plot(tmp2,tmp1(:,3)+yoffset,'r-','LineWidth',6); 
hold off

%% ECDF 2016
for i=1:2
[tmp, ind] = min(abs(time-tson16(1)));
snrall = snr(1:ind); % Only pre-exposure baseline data
dur = diff(find(snrall>dt(i))) - 1;
dur = dur(dur>4); % (>10 mins only)
dur=dur*2.5/60;
    
% Find duration of silence after sonar
nbwDet1 = snr>dt(i);
tdet = time(nbwDet1);
idet = find(nbwDet1==1);
ibefore = find(tdet<tson16(2),1,'last');
sil =(diff(idet(ibefore:ibefore+1))-1)*2.5/60; % duration of observed click-absent period [h]

% Plot ECDFs
subplot(2,2,4); 
[f,x] = ecdf(dur);
hold on; plot(x,f,'-','Color',col(i,:),'LineWidth',3)

% PLot 95%-tile and P of observation
x95 = interp1(f,x,0.95);
hb=plot([x95 x95],[0 0.95],'-','Color',col(i,:),'LineWidth',1);
set(get(get(hb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

end
plot([sil sil],[0 1],'r-','LineWidth',2);
hold off

ylabel('Proportion <= x')
xlabel('Duration of click-absent period (h)')
% legend('DT=5 dB','DT=8 dB','Sonar exposure',4)

