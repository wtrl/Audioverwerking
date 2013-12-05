%% MUSIC_narrowband.m
% group number = 4
% group members: Wouter Lanneer & Philippe de Potter de ten Broeck

clear; clc;
addpath(genpath('..'));
%% Load RIRs and set parameters
load('Computed_RIRs.mat');

% check RIR sampling frequency; must be 44.1kHz
if(fs_RIR ~= 44100)
    error('fs RIR = %d, terwijl het 44.1kHz moet zijn', fs_RIR);
end

% number of microphones
nrOfMics = size(m_pos,1);
% number of sources
N_speech = size(s_pos,1);
N_noise = size(v_pos,1);

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';
speechfilename{2,1} = 'speech2.wav';

% noisefilename = cell(N_noise,3);
% noisefilename{1,1} = 'Babble_noise1.wav';

% length of the recorded microphone signals in seconds
length_recmicsig = 10;
% number of samples in the recording length
nrOfSamples = fs_RIR*length_recmicsig;

%% Create the microphone signals
speechMatrix = zeros(nrOfSamples,N_speech);
noiseMatrix = zeros(nrOfSamples,N_noise);
mic = zeros(nrOfSamples,nrOfMics);
noiseOnMic = zeros(nrOfSamples,nrOfMics);

% read in all speech and noise audiofiles
% resampling of speech and noise files and cut off to number of samples
for i = 1:N_speech
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:nrOfSamples);
end

% for i = 1:N_noise
%     [noisefilename{i,2}, noisefilename{i,3}] = audioread(noisefilename{i,1});
%     noiseTmp = resample(noisefilename{i,2},fs_RIR,noisefilename{i,3});
%     noiseMatrix(:,i) = noiseTmp(1:nrOfSamples);
% end
% noiseMatrix = wgn(nrOfSamples,N_noise,1,'dBm'); % use wgn instead of noise signal

% filter operation
for i = 1:nrOfMics
    % filter speech
    for j = 1:N_speech
        mic(:,i) = mic(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
    % add noise
    %     for j = 1:N_noise
    %         noiseOnMic(:,i) = noiseOnMic(:,i) + fftfilt(RIR_noise(:,i,j),noiseMatrix(:,j));
    %     end
    %     mic(:,i) = mic(:,i) + noiseOnMic(:,i);
end

% calculate noise covariance matrix
% noiseCovarianceMatrix = noiseOnMic'*noiseOnMic;

% save microphone signals and sampling rate in file
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

%% Short Time Fourier Transform
% window length
window = 512;
% 50% overlap
noverlap = window/2;
% number of points for the DFT
nfft = 1024;
% number of time frames
nrOfFrames = fix((length(mic(:,1))-noverlap)/(window-noverlap));
% number of frequency bins
nrOfBins = nfft/2 + 1;

S = zeros(nrOfBins,nrOfMics,nrOfFrames);
P = zeros(nrOfBins,nrOfMics,nrOfFrames);
% STFT
for i = 1:nrOfMics
    x = mic(:,i);
    [S(:,i,:),F,~,P(:,i,:)] = spectrogram(x,window,noverlap,nfft,fs_RIR);
end

<<<<<<< HEAD
% % alternative way of calculating power in each frequency bin
% P_total = zeros(nrOfBins,nrOfMics);
% % average out time
% for i = 1:nrOfMics
%     for j =1:nrOfFrames
%      P_total(:,i) = P_total(:,i) + P(:,i,j);
%     end
% end
% % average out mics
% P_total = mean(P_total,2);
% 
% % S_avg = zeros(NrOfBins, NrOfMics);
% % for t=1:NrOfFrames
% %     S_avg = S_avg + S(:,:,t);
% % end
% 
% 
% [~,maxFreqBin2] = max(P_total);
% maxOmega2 = 2*pi*F(maxFreqBin2);

%% DOA
binCorrMatrix = zeros(nrOfMics, nrOfMics, nrOfBins);
for t = 1:nrOfFrames
    for i = 1:nrOfBins
        binCorrMatrix(:,:,i) = binCorrMatrix(:,:,i) + S(i,:,t)'*S(i,:,t);
    end
end

PSD = zeros(nrOfBins,1);
for i = 1:nrOfBins
    PSD(i) = trace(binCorrMatrix(:,:,i));
=======
%% Select frequency bin with the highest power
P_avg = zeros(nrOfBins,1);
% average out over time and over all the microphones
for i = 1:nrOfMics % average over all the microphones
    P_avg = P_avg + sum(P(:,i,:),3); % the sum command averages over all time frames
>>>>>>> origin/master
end

[maxP_avg,maxFreqBin] = max(P_avg);
maxOmega = 2*pi*F(maxFreqBin);

% plot the PSD
figure;plot(P_avg);
title('Signal power in each frequency bin');xlabel('frequency bins');ylabel('signal power');
hold on; stem(maxFreqBin,maxP_avg,'r');

%% Evaluate the pseudospectrum
% calculate correlation matrix for the frequency bin with highest power
y = zeros(nrOfMics,nrOfFrames);
for i = 1:nrOfMics
    y(i,:) = S(maxFreqBin,i,:);
end
R = y*y';

% % alternative calculation of PSD
% % this however requires to calculate the correlation matrix for each
% % frequency bin
% PSD = zeros(nrOfBins,1);
% for i = 1:nrOfBins
%     PSD(i) = trace(binCorrMatrix(:,:,i));
% end
%
% [maxPSD,maxFreqBin2] = max(abs(PSD));
% maxOmega2 = 2*pi*F(maxFreqBin2);
% figure; plot(abs(PSD));
% title('Signal power in each frequency bin');xlabel('frequency bins');ylabel('signal power');
% hold on; stem(maxFreqBin2,maxPSD,'r');

[V,D] = eigs(R,nrOfMics);
E = V(:,N_speech+1:end);

% distances between microphones in the array relative to first microphone
dm = m_pos(:,2)-m_pos(1,2);

% angle theta; must be in radians
theta = 0:0.5*pi/180:pi; % in radians
g_theta = zeros(nrOfMics, length(theta));
c=340;
for i = 1:length(theta)
    g_theta(:,i) = exp( -1i.*maxOmega*(-dm.*cos(theta(i))./c) );
end

% calculate pseudospectrum for each angle for the maximal frequency bin
Pseudospectrum = zeros(length(theta),1);
for i = 1:length(theta)
    Pseudospectrum(i) =  1/( (g_theta(:,i)')*E*(E')*g_theta(:,i) );
end

%% DOA estimation
% find peaks in the pseudospectrum
[peaks,locs] = findpeaks(abs(Pseudospectrum),'SORTSTR','descend');
maxThetaBins = locs(1:N_speech);
maxsPseudo = peaks(1:N_speech);
DOA_est = (180/pi)*theta(maxThetaBins);
figure; plot(0:0.5:180,abs(Pseudospectrum));

% plot pseudospectrum together with selected peaks
title('Pseudospectrum in function of angle theta');xlabel('theta [degrees]');ylabel('pseudospectrum');
hold on; stem(DOA_est,maxsPseudo,'r');

% save estimate of DOA in file
savefile = 'DOA_est.mat';
save(savefile, 'DOA_est');