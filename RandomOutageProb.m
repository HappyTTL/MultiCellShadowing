%%%%%%%%% Simulation of outage probability given different densities of
%%%%%%%%% Base Station %%%%%%%

%%%%%%%%% Generate a shadow fading field first %%%%%%%%
clear all;
clc;


%%%%%B is a 3000*3000 dimensional matrix represented a correlated shadow
%%%%%fading field
load 'ShadowField9.mat';
ShadowField = B;


%%%%%%%%generate different PPP with different densities
Min_density = 3;
Max_density = 3;
SINR_Threshold = 0; %%%%%SINR threshold is 10db
alpha = 4;
N0 = -104;

sigma = 0;
for Lambda = Min_density:1:Max_density
    %%%%%%%%calculate outage probability
    P0 = 40/Lambda;
    %%%%%%First generate base stations
    Sample_Size = 100000;
    Max_SINR = zeros(Sample_Size,1);
    Max_SIR = zeros(Sample_Size,1);
    Max_SINR_Index = zeros(Sample_Size,1);
    Max_SIR_Index = zeros(Sample_Size,1);
    r = 1000;
    disp('Lambda');
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
     BaseStationPoints = Lambda;
     BaseStationPosition_X = rand(Lambda).*r;
     BaseStationPosition_Y = rand(Lambda).*r;
     BaseStationPosition = [BaseStationPosition_X-r/2 BaseStationPosition_Y-r/2];
        
%             figure(1);
%             plot(BaseStationPosition(:, 1), BaseStationPosition(:, 2), '.');
%             %axis([0 1000 0 1000]);
%             %     hold on;
%             %     plot(BaseStationPoissonProc2(:, 1), BaseStationPoissonProc2(:, 2), '*');
%             %     axis([0 1000 0 1000]);
%             hold on;
%             voronoi(BaseStationPosition(:, 1), BaseStationPosition(:, 2));
        %         filename = ['BSPosition' num2str(n) '.mat'];
        %         save(filename, 'BaseStationPosition');
        
        MU = [0 0];
        N = size(BaseStationPosition);
        Dis = zeros(N(1),1);
        %%%%%%%%%%%calculate pathloss
        PL = zeros(N(1),1);
        PL_DB = zeros(N(1),1);
        for i = 1 : N(1)
            Dis(i) = sqrt(BaseStationPosition(i,1)^2 + BaseStationPosition(i,2)^2);
            PL(i) = Dis(i)^(-alpha);
            PL_DB(i) = 10*log10(PL(i));
        end
        
        
        %%%%%%%%%%calculate shadow fading%%%%%%%%%%%%%
        Shadowing_DB = zeros(N(1),1);
        Shadowing = zeros(N(1),1);
        for i = 1 : N(1)
            Shadowing_DB(i) = B(round(BaseStationPosition(i,1))+1500, round(BaseStationPosition(i,2))+1500)*sigma;
            Shadowing(i) = 10^(Shadowing_DB(i)/10);
        end
        
        %%%%%%%%%%calculate Rayleigh fading %%%%%%%%
        Rayleigh_DB = zeros(N(1),1);
        Rayleigh = zeros(N(1),1);
        for i = 1 : N(1)
            x = sqrt(1/2).*randn(1);
            y = sqrt(1/2).*randn(1);
            Rayleigh(i) = x.^2 + y.^2;
            Rayleigh_DB(i) = 10*log10(Rayleigh(i));
        end
        
        
        RxPower = zeros(N(1),1);
        for i = 1 : N(1)
%             RxPower(i) = 10^(P0/10) * 0.001 * PL(i) * Shadowing(i) * Rayleigh(i);
            RxPower(i) = 10^(P0/10) * 0.001 * PL(i) * Shadowing(i);
        end
        
        SINR = zeros(N(1),1);
        SIR = zeros(N(1),1);
        SINR_DB = zeros(N(1),1);
        SIR_DB = zeros(N(1),1);
        for i = 1 : N(1)
            SINR(i) = RxPower(i)/(sum(RxPower)-RxPower(i)+10^(N0/10));
            SIR(i) = RxPower(i)/(sum(RxPower)-RxPower(i));
            SINR_DB(i) = 10*log10(SINR(i));
            SIR_DB(i) = 10*log10(SIR(i));
        end
        [Max_SINR(n), Max_SINR_Index(n)] = max(SINR_DB);
        [Max_SIR(n), Max_SIR_Index(n)] = max(SIR_DB);
        
    end
    filename = ['Max_SINR' num2str(Lambda) '.mat'];
    save(filename, 'Max_SINR');
    filename = ['Max_SIR' num2str(Lambda) '.mat'];
    save(filename, 'Max_SIR');
    figure(1);
    cdfplot(Max_SINR);
    hold on;
    figure(2);
    cdfplot(Max_SIR);
    hold on;
end
