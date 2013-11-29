%% DOA_corr.m
% group number = 4
% group members: Wouter Lanneer & Philippe de Potter de ten Broeck

clear; clc;
%% Load RIRs and set parameters
load('Computed_RIRs.mat');

% number of microphones
nrOfMics = size(m_pos,1);
% number of sources
N_speech = size(s_pos,1);
N_noise = size(v_pos,1);

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';

% noisefilename = cell(N_noise,3);
% noisefilename{1,1} = 'Babble_noise1.wav';

% length of the recorded microphone signals in seconds
length_recmicsig = 5;
% number of samples in the recording length
nrOfSamples = fs_RIR*length_recmicsig;

%% Calculate ground truth TDOA
indicesMax = zeros(nrOfMics,1);
for i=1:nrOfMics
    [~,indicesMax(i)] = max(RIR_sources(:,i));
end
ground_truth_TDOA = indicesMax(1)-indicesMax(2);
display(ground_truth_TDOA);

%% Create the microphone signals
speechMatrix = zeros(nrOfSamples,N_speech);
noiseMatrix = zeros(nrOfSamples,N_noise);
mic = zeros(nrOfSamples,nrOfMics);

% read in all speech and noise audiofiles
% resampling of speech and noise files and cut off to number of samples
for i = 1:N_speech
    %     speechMatrix = wgn(NrOfSamples,1,1);
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:nrOfSamples);
end

% for i = 1:N_noise
%     [noisefilename{i,2}, noisefilename{i,3}] = audioread(noisefilename{i,1});
%     noiseTmp = resample(noisefilename{i,2},fs_RIR,noisefilename{i,3});
%     noiseMatrix(:,i) = noiseTmp(1:NrOfSamples);
% end

% filter operation
for i = 1:nrOfMics
    % filter speech
    for j = 1:N_speech
        mic(:,i) = mic(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
    % add noise
    %     for j = 1:N_noise
    %         mic(:,i) = mic(:,i) + fftfilt(RIR_noise(:,i,j),noiseMatrix(:,j));
    %     end
end

% save microphone signals and sampling rate in file
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

%% Plot the microphone signals
figure; plot(mic(:,1),'r'); title('The 2 microphone signals');
hold on
plot(mic(:,2),'b');
hold off

%% Cross-correlation-based TDOA estimation 
% select representative segment of the signal and plot it
mic1_crossc = mic(:,1);
figure; plot(mic1_crossc); title('Representative signal segment');
% calculate time-domain cross-correlation function
c = xcorr(mic1_crossc,mic(:,2));
% find the index of the maximum cross-correlation peak
[max_c,estimated_TDOA] = max(c);
estimated_TDOA = estimated_TDOA-nrOfSamples;
display(estimated_TDOA);
% plot cross-correlation function together with the ground truth TDOA
figure; plot(c,'b'); title('Time domain cross-correlation function with ground truth estimate');
hold on; stem(nrOfSamples+ground_truth_TDOA,max_c,'r');
% difference between the estimated TDOA (based on cross-correlation) and
% the ground truth TDOA
diff_TDOA = abs(estimated_TDOA-ground_truth_TDOA);
display(diff_TDOA);

%% Cross-correlation-based DOA estimation 
% speed of sound in air in m/s
c = 340;
% inter-microphone distance in m
d = norm(m_pos(1,:)-m_pos(2,:));
% DOA estimate in degrees between 0-180 degrees
arg = estimated_TDOA*c/(d*fs_RIR);
if(arg > 1)
    DOA_est = 0;
elseif(arg < -1)
    DOA_est = 180;
else
    DOA_est = acos(estimated_TDOA*c/(d*fs_RIR))*180/pi;
end

% save estimate of DOA in file
savefile = 'DOA_est.mat';
save(savefile, 'DOA_est');