%% MUSIC_wideband.m
% group number = 4
% group members: Wouter Lanneer & Philippe de Potter de ten Broeck

clear; clc;
addpath(genpath('..'));
%% Load RIRs and set parameters
load('./Computed_RIRs.mat');

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
% Hann window
L = 1024;
window = hann(L);

% Hamming window
% window = L;

% 50% overlap
noverlap = L/2;
% number of points for the DFT
nfft = L;
% number of time frames
nrOfFrames = fix((length(mic(:,1))-noverlap)/(L-noverlap));
% number of frequency bins
nrOfBins = nfft/2 + 1;

S = zeros(nrOfBins,nrOfMics,nrOfFrames);
P = zeros(nrOfBins,nrOfMics,nrOfFrames);
% STFT
for i = 1:nrOfMics
    x = mic(:,i);
    [S(:,i,:),F,~,P(:,i,:)] = spectrogram(x,window,noverlap,nfft,fs_RIR);
end

%% Evaluate the geometrically averaged pseudospectrum
% calculate correlation matrix for each frequency bin and averaged over all
% time frames
binCorrMatrix = zeros(nrOfMics, nrOfMics, nrOfBins);
y = zeros(nrOfMics,nrOfFrames);
for i = 1:nrOfBins
    for j = 1:nrOfMics
        y(j,:) = S(i,j,:);
    end
    binCorrMatrix(:,:,i) = y*y';
end

% E for each frequency bin
E = zeros(nrOfMics,nrOfMics-N_speech,nrOfBins);
for i = 1:nrOfBins
    [V,D] = eigs(binCorrMatrix(:,:,i),nrOfMics);
    E(:,:,i) = V(:,N_speech+1:end);
end

% distances between microphones in the array relative to first microphone
dm = m_pos(:,2)-m_pos(1,2);

% angle theta; must be in radians
theta = 0:0.5*pi/180:pi; % in radians
omega = 2*pi*F;
g_theta = zeros(nrOfMics,length(theta),nrOfBins);
c=340;
for i = 1:nrOfMics
    for j = 1:nrOfBins
        g_theta(i,:,j) = exp( -1i.*omega(j)*(-dm(i).*cos(theta)./c) );
    end
end

% calculate pseudospectrum for each angle and for each frequency bin
Pseudospectra = zeros(length(theta),nrOfBins);
for i = 1:length(theta)
    for j = 1:nrOfBins
        Pseudospectra(i,j) =  1/( (g_theta(:,i,j)')*E(:,:,j)*(E(:,:,j)')*g_theta(:,i,j) );
    end
end

% geometric average of pseudospectra
logPseudo = (1/(L/2-1))*sum(log(abs(Pseudospectra(:,2:L/2))),2);
Pseudospectrum = exp(logPseudo);

%% DOA estimation
% find peaks in the geometrically averaged pseudospectrum
[peaks,locs] = findpeaks(abs(Pseudospectrum),'SORTSTR','descend');
maxThetaBins = locs(1:N_speech);
maxsPseudo = peaks(1:N_speech);
DOA_est = (180/pi)*theta(maxThetaBins);

% plot geometrically averaged pseudospectra with selected peaks
figure;subplot(212);plot(0:0.5:180,abs(Pseudospectrum));
title('Pseudospectrum in function of angle theta');xlabel('theta [degrees]');ylabel('pseudospectrum');
hold on; stem(DOA_est,maxsPseudo,'r');hold off;
subplot(211);
plot(0:0.5:180,abs(Pseudospectra(:,2:L/2)));
hold on;
hold off;

% save estimate of DOA in file
savefile = 'DOA_est.mat';
save(savefile, 'DOA_est');