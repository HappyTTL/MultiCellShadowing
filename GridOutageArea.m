clear all;
clc


SINR_Threshold = 0; %%%%%SINR threshold is 10db
alpha = 4;
N0 = -104;
P0 = 30;
sigma = 8;
%%%%%
%% Generate a Random Layout of density Lambda
Lambda = 1;
r = 1000;
n = sqrt(Lambda);
BaseStationPosition = GridBaseLayout(r,n);
% BaseStationPoints = Lambda;
% BaseStationPosition_X = rand(Lambda) .* r;
% BaseStationPosition_Y = rand(Lambda) .* r;
% BaseStationPosition = [BaseStationPosition_X - r/2 BaseStationPosition_Y - r/2];
figure(1);
plot(BaseStationPosition(:, 1), BaseStationPosition(:, 2), '.');
axis([-500 500 -500 500]);
% hold on;
% voronoi(BaseStationPosition(:, 1), BaseStationPosition(:, 2));
% title('Random Layout');
%% Load a correlated shadow fading field and calculate the path loss of 
%%% with shadow fading and without shadow fading

SF = zeros(Lambda,3000,3000);
load('ShadowField9.mat');
SF(1,:,:) = B;
load('ShadowField8.mat');
SF(2,:,:) = B;
load('ShadowField7.mat');
SF(3,:,:) = B;
load('ShadowField6.mat');
SF(4,:,:) = B;
load('ShadowField5.mat');
SF(5,:,:) = B;
load('ShadowField4.mat');
SF(6,:,:) = B;
load('ShadowField3.mat');
SF(7,:,:) = B;
load('ShadowField2.mat');
SF(8,:,:) = B;
load('ShadowField1.mat');
SF(9,:,:) = B;


Max_SINR_NoShadow_DB = zeros(r,r);
Max_SINR_WithShadow_DB = zeros(r,r);
Max_SINR_IIDShadow_DB = zeros(r,r);
SINR_Marker_NoShadow = (zeros(r,r) + 1) .* 255;
SINR_Marker_WithShadow = (zeros(r,r) + 1) .* 255;
SINR_Marker_IIDShadow = (zeros(r,r) + 1) .* 255;
SINR_Marker_NoShadow_2 = (zeros(r,r) + 1) .* 255;
SINR_Marker_WithShadow_2 = (zeros(r,r) + 1) .* 255;
Index_NoShadow = zeros(r,r);
Index_WithShadow = zeros(r,r);
Index_IIDShadow = zeros(r,r);
%% Calculate pathloss for each point in a r*r area
for m = 1 : r
    for n = 1 : r
        Dis = zeros(Lambda,1);
        PL = zeros(Lambda,1);
        PL_DB = zeros(Lambda,1);
        for i = 1 : Lambda
            Dis(i) = sqrt((BaseStationPosition(i,1) - m + r/2)^2 + (BaseStationPosition(i,2) - n + r/2)^2);
            PL(i) = Dis(i)^(-alpha);
            PL_DB(i) = 10 * log10(PL(i));
        end
        
        %%%%%Calculate relative position from point (m,n) to every base
        %%%%%station
        
        Relative_Pos = zeros(Lambda,2);
        ShadowFading = zeros(Lambda,1);
        IID_ShadowFading = zeros(Lambda,1);
        ShadowFading_DB = zeros(Lambda,1);
        IID_ShadowFading_DB = zeros(Lambda,1);
        for i = 1 : Lambda
            Relative_Pos(i,:) = [round(m - r/2 - BaseStationPosition(i,1))  round(n - r/2 - BaseStationPosition(i,2))];
            ShadowFading_DB(i) = SF(i, Relative_Pos(i,1) + 1500, Relative_Pos(i,2) + 1500) * sigma;
            ShadowFading(i) = 10^(ShadowFading_DB(i)/10);
            IID_ShadowFading_DB(i) = normrnd(0, sigma);
            IID_ShadowFading = 10^(IID_ShadowFading_DB(i)/10);
        end
        
        %%%%%%%%%%calculate Rayleigh fading %%%%%%%%
        Rayleigh_DB = zeros(Lambda,1);
        Rayleigh = zeros(Lambda,1);
        for i = 1 : Lambda
            x = sqrt(1/2).*randn(1);
            y = sqrt(1/2).*randn(1);
            Rayleigh(i) = x.^2 + y.^2;
            Rayleigh_DB(i) = 10*log10(Rayleigh(i));
        end
        
        %%%%%Calculate Received Power
        RxPower_WithShadow = zeros(Lambda,1);
        RxPower_NoShadow = zeros(Lambda,1);
        RxPower_IIDShadow = zeros(Lambda,1);
        for i = 1 : Lambda
            RxPower_WithShadow(i) = 10^(P0/10) * 0.001 * PL(i) * ShadowFading(i) * Rayleigh(i);
            RxPower_NoShadow(i) = 10^(P0/10) * 0.001 * PL(i) * Rayleigh(i);
            RxPower_IIDShadow(i) = 10^(P0/10) * 0.001 * PL(i) * IID_ShadowFading(i) * Rayleigh(i);
        end
        
        SINR_NoShadow = zeros(Lambda,1);
        SINR_WithShadow = zeros(Lambda,1);
        SINR_NoShadow_DB = zeros(Lambda,1);
        SINR_WithShadow_DB = zeros(Lambda,1);
        SINR_IIDShadow = zeros(Lambda,1);
        SINR_IIDShadow_DB = zeros(Lambda,1);
        for i = 1 : Lambda
            SINR_NoShadow(i) = RxPower_NoShadow(i)/(sum(RxPower_NoShadow)-RxPower_NoShadow(i)+10^(N0/10));
            SINR_WithShadow(i) = RxPower_WithShadow(i)/(sum(RxPower_WithShadow)-RxPower_WithShadow(i)+10^(N0/10));
            SINR_IIDShadow(i) = RxPower_IIDShadow(i)/(sum(RxPower_IIDShadow)-RxPower_IIDShadow(i)+10^(N0/10));
            SINR_IIDShadow_DB(i) = 10*log10(SINR_IIDShadow(i));
            SINR_NoShadow_DB(i) = 10*log10(SINR_NoShadow(i));
            SINR_WithShadow_DB(i) = 10*log10(SINR_WithShadow(i));
        end
        

        [Max_SINR_NoShadow_DB(m,n), Index_NoShadow(m,n)] =  max(SINR_NoShadow_DB);
        [Max_SINR_WithShadow_DB(m,n), Index_WithShadow(m,n)] =  max(SINR_WithShadow_DB);
        [Max_SINR_IIDShadow_DB(m,n), Index_IIDShadow(m,n)] =  max(SINR_IIDShadow_DB);
        
        if Max_SINR_NoShadow_DB(m,n) < SINR_Threshold
            SINR_Marker_NoShadow(m,n) = 0;
        end
        
        if Max_SINR_WithShadow_DB(m,n) < SINR_Threshold
            SINR_Marker_WithShadow(m,n) = 0;
        end
        
        if Max_SINR_IIDShadow_DB(m,n) < SINR_Threshold
            SINR_Marker_IIDShadow(m,n) = 0;
        end
        
        if Max_SINR_NoShadow_DB(m,n) < SINR_Threshold - 5
            SINR_Marker_NoShadow_2(m,n) = 0;
        end
        
        if Max_SINR_WithShadow_DB(m,n) < SINR_Threshold - 5
            SINR_Marker_WithShadow_2(m,n) = 0;
        end
        
    end
end

figure(2);
imagesc(SINR_Marker_IIDShadow);
axis([0 r 0 r]);
colormap(gray);
title('Outage Area with IID Shadow Fading');
figure(3);
subplot(1,2,1);
imagesc(SINR_Marker_WithShadow);
axis([0 r 0 r]);
colormap(gray);
title('Outage Area with Shadow Fading'); 
subplot(1,2,2);
imagesc(SINR_Marker_IIDShadow);
axis([0 r 0 r]);
colormap(gray);
title('Outage Area with IID Shadow Fading'); 
figure(4);
imagesc(Max_SINR_NoShadow_DB);
colormap(jet);
title('Max SINR without Shadow Fading');
figure(5);
imagesc(Max_SINR_WithShadow_DB);
colormap(jet);
title('Max SINR with Shadow Fading');


        
        
        
        
        
        
        