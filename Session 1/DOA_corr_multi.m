clear;

load('Computed_RIRs.mat');

%Number of mics
NrOfMics = size(m_pos,1);

%Number of files
N_speech = size(s_pos,1);

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
for i=1:N_speech
    speechfilename{i,1} = 'speech1.wav';
end

% length of the recorded microphone signals in seconds
length_recmicsig = 5;

NrOfSamples = fs_RIR*length_recmicsig;

%Ground truth TDOA
ground_truth_TDOA = zeros(N_speech,1);

indicesMax = zeros(NrOfMics,N_speech);
for j = 1:N_speech
    for i = 1:NrOfMics
        [~,indicesMax(i,j)] = max(RIR_sources(:,i,j));
    end
end
for i = 1:N_speech
    ground_truth_TDOA(i) = indicesMax(1,i)-indicesMax(2,i);
end

speechMatrix = zeros(NrOfSamples,N_speech);
mic = zeros(NrOfSamples,NrOfMics,N_speech);

% read in all speech audiofiles
for i = 1:N_speech
    %     speechMatrix = wgn(NrOfSamples,1,1);
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:NrOfSamples);
end

%filteroperation
% creating the microphone signals when only one source is active at a time
for j = 1:N_speech
    for i = 1:NrOfMics
        %filter speech
        mic(:,i,j) = fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
end

%save files;
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

% figure; plot(mic(:,1),'r'); title('The 2 microphone signals');
% hold on
% plot(mic(:,2),'b');
% hold off

% cross correlation method for all speech sources
estimated_TDOA = zeros(N_speech,1);

for i = 1:N_speech
    c = xcorr(mic(:,1,i),mic(:,2,i));
    [max_c,estimated_TDOA(i)] = max(c);
    estimated_TDOA(i) = estimated_TDOA(i)-NrOfSamples;
end

% figure; plot(c,'b'); title('Time domain cross-correlation function with ground truth estimate');
% hold on; stem(NrOfSamples+ground_truth_TDOA,max_c,'r');

diff_TDOA = (estimated_TDOA-ground_truth_TDOA)

%%%%%%%%%
%   DOA %
%%%%%%%%%

% speed of sound in air in m/s
c = 340;

% inter-microphone distance in m
d = norm(m_pos(1,:)-m_pos(2,:));

% DOA estimate in degrees between 0-180 degrees
arg = estimated_TDOA.*(c/(d*fs_RIR));
arg(arg > 1) = 1;
arg(arg < -1) = -1;
DOA_est = acos(arg)*180/pi;

%save estimate of DOA in file
savefile = 'DOA_est.mat';
save(savefile, 'DOA_est');

% figure;
% t = -1:0.01:1;
% plot(t,acos(t)); title('Acos function');