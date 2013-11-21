clear;

load('Computed_RIRs.mat');






figure;
RIR = RIR_sources(:,1);
plot(RIR, 'b');
hold on
plot(RIR_sources(:,2),'r');
title('RIR from source 2 to microphone 3');