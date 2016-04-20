clear all;
clc

%% load data
load('Max_SINR_Outage9.mat');

Max_SINR(find(Max_SINR >= -10)) = 1;
Max_SINR(find(Max_SINR < -10)) = 0;

S = size(Max_SINR);
NewMax_SINR = [Max_SINR'; zeros(1,S(1)) + 1]';
OneDim_Max_SINR = reshape(NewMax_SINR',[],1);
LengthOneSeqCount = 0;
S_All = size(OneDim_Max_SINR);
for i = 2 : S_All(1)-1
    if OneDim_Max_SINR(i-1) == 1 && OneDim_Max_SINR(i) == 0 && OneDim_Max_SINR(i+1) == 1
        LengthOneSeqCount = LengthOneSeqCount + 1; 
    end
end
LengthOneSeq = zeros(LengthOneSeqCount,1) + 1;
Seq = findseq(Max_SINR,2);

S = size(Seq);

Outage_Seq_All = zeros(S(1),1);
j = 1;
for i = 1 : S(1)
    if Seq(i,1) == 0
        Outage_Seq_All(j) = Seq(i,4);
        j = j + 1;
    end
end

Outage_Seq = [LengthOneSeq; Outage_Seq_All(find(Outage_Seq_All ~= 0))];

figure(1);
cdfplot(Outage_Seq);
hold on;
        