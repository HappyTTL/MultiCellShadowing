%%%%%%%%% Simulation of outage probability given different densities of
%%%%%%%%% Base Station %%%%%%%

%%%%Author: Tingting Lu
%%%%Date: 3/21/2016

%%%%%%BS and MU positions%%%%%%
%%Two cell size: 1000m and 500m

%% 
clear all
clc
bestChannel = 1;
nearestBaseStation = 0;
r = 1000;
n = 3;
Location = GridBaseLayout(r,n);
figure(1);
scatter(Location(:,1), Location(:,2));
axis([-r/2 r/2 -r/2 r/2]);
voronoi(Location(:, 1), Location(:, 2));
xlim([-r/2 r/2]);
ylim([-r/2 r/2]);
title('Grid Layout');
%%
SNRThreshold = 0;
SINR_Low = -30;
SINR_High = 30;
N = 100000;
SINR_Mark_Corr = zeros(SINR_High-SINR_Low+1,N);
SINR_Mark_No = zeros(SINR_High-SINR_Low+1,N);
SINR_Mark_IID = zeros(SINR_High-SINR_Low+1,N);
SIR_Mark_Corr = zeros(SINR_High-SINR_Low+1,N);
SIR_Mark_No = zeros(SINR_High-SINR_Low+1,N);
SIR_Mark_IID = zeros(SINR_High-SINR_Low+1,N);
sigma = 8;
alpha = 4;
N0 = -104;
% P0 = N0 + SNRThreshold;
P0 = 30;
%% 


MU_x = zeros(N,1);
MU_y = zeros(N,1);
Dis = zeros(N,n^2);
PL_DB = zeros(N,n^2);
PL = zeros(N,n^2);
S = zeros(N,n^2);
CG_Corr = zeros(N,n^2);
CG_No = zeros(N,n^2);
CG_IID = zeros(N,n^2);
SINR_Corr = zeros(N,1);
SINR_IID = zeros(N,1);
SINR_No = zeros(N,1);
SIR_Corr = zeros(N,1);
SIR_No = zeros(N,1);
SIR_IID = zeros(N,1);
%% 

for i = 1 : N
    if mod(n,2) == 1
        MU_x(i) = rand*r/n - r/2/n;
        MU_y(i) = rand*r/n - r/2/n;
    else
        MU_x(i) = rand*r/n;
        MU_y(i) = rand*r/n;
    end
% 
% figure(2);
% scatter(MU_x, MU_y);

%%%%%First Compute Distance%%%%%%%%

    
    %%%%%%%%%%%%path loss%%%%%%%%%%%%%%%
    for j = 1 : n^2
        Dis(i,j) = sqrt((MU_x(i)-Location(j,1))^2 + (MU_y(i)-Location(j,2))^2);
        %PL_DB1(i,j) = 15.3+37.6*log10(Dis1(i,j));
        PL(i,j) = Dis(i,j)^(-alpha);
        PL_DB(i,j) = 10*log10(PL(i,j));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%Rayleigh Fading%%%%%%%%%%%%%%%%%
    x = sqrt(1/2).*randn(1);
    y = sqrt(1/2).*randn(1);
    hsquare = x.^2 + y.^2;
    
    %%%%%%%%%%%%%%%%%%%%%Shadow Fading%%%%%%%%%%%%%%%%%%%%%
    R = 6;
    Dis_R = zeros(n^2,n^2);
    R0 = zeros(n^2,n^2);
    H_R0 = zeros(n^2,n^2);
    theta = zeros(n^2,n^2);
    H_theta = zeros(n^2,n^2);
    H = zeros(n^2,n^2);
    
    for j = 1:n^2
        for p = j:n^2
            Dis_R(j,p) = sqrt((Location(j,1)-Location(p,1))^2 + (Location(j,2)-Location(p,2))^2);
            R0(j,p) = abs(10*log10(Dis(i,j)/Dis(i,p)));
            H_R0(j,p) = max(0, 1-R0(j,p)/R);
            theta(j,p) = acos((Dis(i,j)^2+Dis(i,p)^2-Dis_R(j,p)^2)/2/Dis(i,j)/Dis(i,p));            %theta2(j,p) = acos((Dis2(i,j)^2+Dis2(i,p)^2-Dis_R2(j,p)^2)/2/Dis2(i,j)/Dis2(i,p));
            if theta(j,p) > pi/3
                H_theta(j,p) = 0;
            else H_theta(j,p) = 1 - theta(j,p)/(pi/3);
            end
            H(j,p) = H_R0(j,p)*H_theta(j,p);
            H(p,j) = H(j,p);
        end
    end
    C = chol(real(H));
    Z = normrnd(0,sigma,[n^2,1]);    %S1(i,:) = C1'*Z;    %S2(i,:) = C2'*Z;
    S(i,:) = 10.^(C*Z./10);
    
    
    CG_Corr(i,:) = P0 * PL(i,:) .* S(i,:) * hsquare;
    CG_IID(i,:) = P0 * PL(i,:) .* 10.^(Z/10)' * hsquare;
    CG_No(i,:) = P0* PL(i,:) * hsquare;
    
    %%%%%%%%%%Best Channel%%%%%%%%%%%%
    if bestChannel == 1
        SIR_Corr(i) = 10*log10(max(CG_Corr(i,:))/(sum(CG_Corr(i,:))-max(CG_Corr(i,:))));
        SINR_Corr(i) = 10*log10(max(CG_Corr(i,:))/(sum(CG_Corr(i,:))-max(CG_Corr(i,:))+10^(N0/10)));
        SIR_No(i) = 10*log10(max(CG_No(i,:))/(sum(CG_No(i,:))-max(CG_No(i,:))));
        SINR_No(i) = 10*log10(max(CG_No(i,:))/(sum(CG_No(i,:))-max(CG_No(i,:))+10^(N0/10)));
        SIR_IID(i) = 10*log10(max(CG_IID(i,:))/(sum(CG_IID(i,:))-max(CG_IID(i,:))));
        SINR_IID(i) = 10*log10(max(CG_IID(i,:))/(sum(CG_IID(i,:))-max(CG_IID(i,:))+10^(N0/10)));
    end
    
    %%%%%%%%%%%%Nearest Base Station%%%%%%%%%%
    if nearestBaseStation == 1
        if mod(n,2) == 1
            SIR_Corr(i) = 10*log10(CG_Corr(i,ceil(n^2/2))/(sum(CG_Corr(i,:))-CG_Corr(i,ceil(n^2/2))));
            SINR_Corr(i) = 10*log10(CG_Corr(i,ceil(n^2/2))/(sum(CG_Corr(i,:))-CG_Corr(i,ceil(n^2/2))+10^(N0/10)));
            SIR_No(i) = 10*log10(CG_No(i,ceil(n^2/2))/(sum(CG_No(i,:))-CG_No(i,ceil(n^2/2))));
            SINR_No(i) = 10*log10(CG_No(i,ceil(n^2/2))/(sum(CG_No(i,:))-CG_No(i,ceil(n^2/2))+10^(N0/10)));
            SIR_IID(i) = 10*log10(CG_IID(i,ceil(n^2/2))/(sum(CG_IID(i,:))-CG_IID(i,ceil(n^2/2))));
            SINR_IID(i) = 10*log10(CG_IID(i,ceil(n^2/2))/(sum(CG_IID(i,:))-CG_IID(i,ceil(n^2/2))+10^(N0/10)));
        else
            SIR_Corr(i) = 10*log10(CG_Corr(i,ceil((n^2+n+1)/2))/(sum(CG_Corr(i,:))-CG_Corr(i,ceil((n^2+n+1)/2))));
            SINR_Corr(i) = 10*log10(CG_Corr(i,ceil((n^2+n+1)/2))/(sum(CG_Corr(i,:))-CG_Corr(i,ceil((n^2+n+1)/2))+10^(N0/10)));
            SIR_No(i) = 10*log10(CG_No(i,ceil((n^2+n+1)/2))/(sum(CG_No(i,:))-CG_No(i,ceil((n^2+n+1)/2))));
            SINR_No(i) = 10*log10(CG_No(i,ceil((n^2+n+1)/2))/(sum(CG_No(i,:))-CG_No(i,ceil((n^2+n+1)/2))+10^(N0/10)));
            SIR_IID(i) = 10*log10(CG_IID(i,ceil((n^2+n+1)/2))/(sum(CG_IID(i,:))-CG_IID(i,ceil((n^2+n+1)/2))));
            SINR_IID(i) = 10*log10(CG_IID(i,ceil((n^2+n+1)/2))/(sum(CG_IID(i,:))-CG_IID(i,ceil((n^2+n+1)/2))+10^(N0/10)));
        end
    end
  
    for m = SINR_Low:1:SINR_High
        
        if SINR_Corr(i) > m
            SINR_Mark_Corr(m - SINR_Low +1, i) = 1;
        end
        if SINR_No(i) > m
            SINR_Mark_No(m - SINR_Low +1, i) = 1;
        end
        if SINR_IID(i) > m
            SINR_Mark_IID(m - SINR_Low +1,i) = 1;
        end
        
        if SIR_Corr(i) > m
            SIR_Mark_Corr(m - SINR_Low +1, i) = 1;
        end
        if SIR_No(i) > m
            SIR_Mark_No(m - SINR_Low +1, i) = 1;
        end
        if SIR_IID(i) > m
            SIR_Mark_IID(m - SINR_Low +1,i) = 1;
        end
    end
end
SINR_Cov_Prob_Corr = zeros(1,SINR_High - SINR_Low +1);
SINR_Cov_Prob_No = zeros(1,SINR_High - SINR_Low +1);
SINR_Cov_Prob_IID = zeros(1,SINR_High - SINR_Low +1);

SIR_Cov_Prob_Corr = zeros(1,SINR_High - SINR_Low +1);
SIR_Cov_Prob_No = zeros(1,SINR_High - SINR_Low +1);
SIR_Cov_Prob_IID = zeros(1,SINR_High - SINR_Low +1);
for m = SINR_Low:1:SINR_High
    
    SIR_Cov_Prob_Corr(m - SINR_Low +1) = sum(SIR_Mark_Corr(m - SINR_Low +1,:));
    SIR_Cov_Prob_No(m - SINR_Low +1) = sum(SIR_Mark_No(m - SINR_Low +1,:));
    SIR_Cov_Prob_IID(m - SINR_Low +1) = sum(SIR_Mark_IID(m - SINR_Low +1,:));
    SINR_Cov_Prob_Corr(m - SINR_Low +1) = sum(SINR_Mark_Corr(m - SINR_Low +1,:));
    SINR_Cov_Prob_No(m - SINR_Low +1) = sum(SINR_Mark_No(m - SINR_Low +1,:));
    SINR_Cov_Prob_IID(m - SINR_Low +1) = sum(SINR_Mark_IID(m - SINR_Low +1,:));

end

save('Grid_SINR_Cov_Prob_25_BestC.mat', 'SINR_Cov_Prob_No','SINR_Cov_Prob_IID','SINR_Cov_Prob_Corr');
save('Grid_SIR_Cov_Prob_25_BestC.mat', 'SIR_Cov_Prob_No','SIR_Cov_Prob_IID','SIR_Cov_Prob_Corr');
figure(2);
% plot((SINR_Low:1:SINR_High), SINR_Cov_Prob_Corr/N, 'r','LineWidth',2);
% hold on;
% plot((SINR_Low:1:SINR_High), SINR_Cov_Prob_No/N, 'g', 'LineWidth',2);
% hold on;
% plot((SINR_Low:1:SINR_High), SINR_Cov_Prob_IID/N, 'k', 'LineWidth',2);
% hold on;
plot((SINR_Low:1:SINR_High), SIR_Cov_Prob_Corr/N, '--r','LineWidth',2);
hold on;
plot((SINR_Low:1:SINR_High), SIR_Cov_Prob_No/N, '--g','LineWidth',2);
hold on;
plot((SINR_Low:1:SINR_High), SIR_Cov_Prob_IID/N, '--k','LineWidth',2);
legend('SIRCovProbCorr','SIRCovProbNo','SIRCovProbIID');
title('Coverage Probability of Grid Layout (Best C 25)');
xlabel('SIR Threshold');
ylabel('Coverage Probability');
axis([SINR_Low SINR_High 0 1]);
grid;
% figure(2);
% cdfplot(real(SINR0));
% hold on;
% cdfplot(real(SINR1));
% hold on;
% cdfplot(real(SINR2));
% hold on;
% cdfplot(real(SINR3));



%%%%Throughput%%%%%%