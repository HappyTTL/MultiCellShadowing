%%%%%%%%% Simulation of outage probability given different densities of
%%%%%%%%% Base Station %%%%%%%

%%%%Author: Tingting Lu
%%%%Date: 3/21/2016

%%%%%%BS and MU positions%%%%%%
%%Two cell size: 1000m and 500m

%% 
clear all
clc

r = 1000;
n = 3;
Location = GridBaseLayout(r,n);
% figure(1);
% scatter(Location(:,1), Location(:,2));
% axis([-r/2 r/2 -r/2 r/2]);
% % voronoi(Location(:, 1), Location(:, 2));
% xlim([-r/2 r/2]);
% ylim([-r/2 r/2]);
% title('Grid Layout');
%%
SNRThreshold = 0;
Sample_Size = 1000000;
sigma = 1;
alpha = 4;
N0 = -104;
% P0 = N0 + SNRThreshold;
P0 = 40;
%% 


MU_x = zeros(Sample_Size,1);
MU_y = zeros(Sample_Size,1);
Dis = zeros(Sample_Size,n^2);
Relative_Dis = zeros(Sample_Size,n^2,2);
PL_DB = zeros(Sample_Size,n^2);
PL = zeros(Sample_Size,n^2);
Shadowing = zeros(Sample_Size,n^2);
Shadowing_DB = zeros(Sample_Size,n^2);
Rayleigh_DB = zeros(Sample_Size,n^2);
Rayleigh = zeros(Sample_Size,n^2);
RxPower = zeros(Sample_Size, n^2);
SINR = zeros(Sample_Size, n^2);
SIR = zeros(Sample_Size, n^2);
SINR_DB = zeros(Sample_Size, n^2);
SIR_DB = zeros(Sample_Size, n^2);
Max_SINR = zeros(Sample_Size,1);
Max_SIR = zeros(Sample_Size,1);
Max_SINR_Index = zeros(Sample_Size,1);
Max_SIR_Index = zeros(Sample_Size,1);
%% 

for i = 1 : Sample_Size
    MU_x(i) = rand * r - r/2;
    MU_y(i) = rand * r - r/2;
end
% 
% figure(2);
% scatter(MU_x, MU_y);

%%%%%First Compute Distance%%%%%%%%

%% load shadowing field
load 'ShadowField1.mat';
ShadowField = zeros(n^2, size(B,1), size(B,2));
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


%% Compute pathloss and shadow fading and Rayleigh fading
    %%%%%%%%%%%%path loss%%%%%%%%%%%%%%%
for i = 1 : Sample_Size
    for j = 1 : n^2
        Dis(i,j) = sqrt((MU_x(i)-Location(j,1))^2 + (MU_y(i)-Location(j,2))^2);
        %PL_DB1(i,j) = 15.3+37.6*log10(Dis1(i,j));
        PL(i,j) = Dis(i,j)^(-alpha);
        PL_DB(i,j) = 10*log10(PL(i,j));
        Relative_Dis(i,j,:) = [MU_x(i) - Location(j,1) MU_y(i) - Location(j,2)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%Rayleigh Fading%%%%%%%%%%%%%%%%%
    for j = 1 : n^2
        x = sqrt(1/2).*randn(1);
        y = sqrt(1/2).*randn(1);
        hsquare = x.^2 + y.^2;
        Rayleigh(i,j) = hsquare;
        Rayleigh_DB(i,j) = 10 * log10(Rayleigh(i,j));
    end
    
    %%%%%%%%%%%%%%%%%%%%%Shadow Fading%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1 : n^2
        Shadowing_DB(i,j) = ShadowField(j, round(Relative_Dis(i,j,1)) + 1500, round(Relative_Dis(i,j,2))+1500) * sigma;
        Shadowing(i,j) = 10^(Shadowing_DB(i,j)/10);
    end
    
    
    %%%%%Receive Power
    
    for j = 1 : n^2
%         RxPower(i,j) = 10^(P0/10) * 0.001 * PL(i,j) * Shadowing(i,j) * Rayleigh(i,j);
        RxPower(i,j) = 10^(P0/10) * 0.001 * PL(i,j) * Shadowing(i,j);
    end
    
    for j = 1 : n^2
        SINR(i,j) = RxPower(i,j)/(sum(RxPower(i,:))-RxPower(i,j)+10^(N0/10));
        SIR(i,j) = RxPower(i,j)/(sum(RxPower(i,:))-RxPower(i,j));
        SINR_DB(i,j) = 10*log10(SINR(i,j));
        SIR_DB(i,j) = 10*log10(SIR(i,j));
    end
    
    [Max_SINR(i), Max_SINR_Index(i)] = max(SINR_DB(i,:));
    [Max_SIR(i), Max_SIR_Index(i)] = max(SIR_DB(i,:));
end

filename = ['Max_SINR_Outage_Grid' num2str(n^2) '.mat'];
save(filename, 'Max_SINR');
filename = ['Max_SIR_Outage_Grid' num2str(n^2) '.mat'];
save(filename, 'Max_SIR');
    
figure(1);
cdfplot(Max_SINR);
hold on;