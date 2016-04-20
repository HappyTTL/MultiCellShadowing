clear all;
clc

Density = 9;
load 'ShadowField1.mat';
ShadowField = zeros(Density, size(B,1), size(B,2));
ShadowField(1,:,:) = B;
load 'ShadowField2.mat';
ShadowField(2,:,:) = B;
load 'ShadowField3.mat';
ShadowField(3,:,:) = B;
load 'ShadowField4.mat';
ShadowField(4,:,:) = B;
load 'ShadowField5.mat';
ShadowField(5,:,:) = B;
load 'ShadowField6.mat';
ShadowField(6,:,:) = B;
load 'ShadowField7.mat';
ShadowField(7,:,:) = B;
load 'ShadowField8.mat';
ShadowField(8,:,:) = B;
load 'ShadowField9.mat';
ShadowField(9,:,:) = B;

%%%%%%%%generate different PPP with different densities
Min_density = 4;
Max_density = 4;
SINR_Threshold = 0; %%%%%SINR threshold is 10db
alpha = 4;
N0 = -104;
P0 = 30;
Sigma = 8;
for Lambda = Min_density:1:Max_density
    %%%%%%%%calculate outage probability for the x axis%%%%%
    
    %%%%%%First generate base stations
    
    r = 1000;
    Sample_Size = round(r/4);
    disp('Lambda');
    
    Max_SINR = zeros(Sample_Size,r);
    Max_SIR = zeros(Sample_Size,r);
    Max_SINR_Index = zeros(Sample_Size,r);
    Max_SIR_Index = zeros(Sample_Size,r);
    for n = 1:Sample_Size
        %%%%%%Generate to PPP base station Layout%%%%%%%

%         BaseStationPoints = poissrnd(Lambda);
%         if BaseStationPoints == 0
%             disp('No BaseStations');
%             break;
%         end
%             
%         %BaseStationPoints2 = poissrnd(BaseStationLambda2);
%         BaseStationPoissonProc = rand(BaseStationPoints, 2).*r;
%         BaseStationPosition = BaseStationPoissonProc-500;
     %%%%%Generate random Base Stations%%%%%%%%%%%
%      BaseStationPoints = Lambda;
%      BaseStationPosition_X = rand(Lambda).*r;
%      BaseStationPosition_Y = rand(Lambda).*r;
%      BaseStationPosition = [BaseStationPosition_X-r/2 BaseStationPosition_Y-r/2];
     Num = sqrt(Lambda);
     BaseStationPosition = GridBaseLayout(r,Num);
     %     figure(1);
        %     plot(BaseStationPosition(:, 1), BaseStationPosition(:, 2), 'o');
        %     %axis([0 1000 0 1000]);
        %     %     hold on;
        %     %     plot(BaseStationPoissonProc2(:, 1), BaseStationPoissonProc2(:, 2), '*');
        %     %     axis([0 1000 0 1000]);
        %     hold on;
        %     voronoi(BaseStationPosition(:, 1), BaseStationPosition(:, 2));
        %         filename = ['BSPosition' num2str(n) '.mat'];
        %         save(filename, 'BaseStationPosition');
        
        MU_Position_X = (-r/2:r/2-1);
        MU_Position_Y = zeros(r,1)-round(r/4)+i;
        MU_Position = [MU_Position_X' MU_Position_Y];
        N = size(BaseStationPosition);
        Dis = zeros(N(1),r);
        %%%%%%%%%%%calculate pathloss
        PL = zeros(N(1),r);
        PL_DB = zeros(N(1),r);
        for i = 1 : N(1)
            for j = 1 : r
                Dis(i,j) = sqrt((BaseStationPosition(i,1)-MU_Position(j,1))^2 + BaseStationPosition(i,2)^2);
                PL(i,j) = Dis(i,j)^(-alpha);
                PL_DB(i,j) = 10*log10(PL(i,j));
            end
        end
        disp('pathloss');
        
        %%%%%%%%%%calculate shadow fading%%%%%%%%%%%%%
        Shadowing_DB = zeros(N(1),r);
        Shadowing = zeros(N(1),r);
        for i = 1 : N(1)
            for j = 1 : r
                Shadowing_DB(i,j) = ShadowField(i,round(BaseStationPosition(i,1))+1500, round(BaseStationPosition(i,2))+1500)*Sigma;
                Shadowing(i,j) = 10^(Shadowing_DB(i,j)/10);
            end
        end
        disp('shadowfading')
        %%%%%%%%%%calculate Rayleigh fading %%%%%%%%
        Rayleigh_DB = zeros(N(1),r);
        Rayleigh = zeros(N(1),r);
        for i = 1 : N(1)
            for j = 1 : r
                x = sqrt(1/2).*randn(1);
                y = sqrt(1/2).*randn(1);
                Rayleigh(i,j) = x.^2 + y.^2;
                Rayleigh_DB(i,j) = 10*log10(Rayleigh(i,j));
            end
        end
        
        
        RxPower = zeros(N(1),r);
        for i = 1 : N(1)
            for j = 1 : r
                RxPower(i,j) = 10^(P0/10) * 0.001 * PL(i,j) * Shadowing(i,j) * Rayleigh(i,j);
            end
        end
        
        SINR = zeros(N(1),r);
        SIR = zeros(N(1),r);
        SINR_DB = zeros(N(1),r);
        SIR_DB = zeros(N(1),r);
        for i = 1 : N(1)
            for j = 1 : r
                SINR(i,j) = RxPower(i,j)/(sum(RxPower(:,j))-RxPower(i,j)+10^(N0/10));
                SIR(i,j) = RxPower(i,j)/(sum(RxPower(:,j))-RxPower(i,j));
                SINR_DB(i,j) = 10*log10(SINR(i,j));
                SIR_DB(i,j) = 10*log10(SIR(i,j));
            end
        end
        for j = 1 : r
            [Max_SINR(n,j), Max_SINR_Index(n,j)] = max(SINR_DB(:,j));
            [Max_SIR(n,j), Max_SIR_Index(n,j)] = max(SIR_DB(:,j));
        end
        
    end
    filename = ['Max_SINR_Outage_Grid' num2str(Lambda) '.mat'];
    save(filename, 'Max_SINR');
    filename = ['Max_SIR_Outage_Grid' num2str(Lambda) '.mat'];
    save(filename, 'Max_SIR');
%     figure(1);
%     plot((1:j), Max_SINR);
%     hold on;
%     figure(2);
%     plot((1:j), Max_SIR);
%     hold on;
end
