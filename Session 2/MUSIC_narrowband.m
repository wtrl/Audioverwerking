clear; clc;
addpath(genpath('..'));

%% Compute all microphone signals

% user defined file names of the WAV files to use for sources and noise
load('Computed_RIRs.mat');

if(fs_RIR ~= 44100)
    error('fs RIR = %d, terwijl het 44.1kHz moet zijn', fs_RIR);
end

% Number of files
N_speech = size(s_pos,1);
N_noise = size(v_pos,1);

speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';
%speechfilename{2} = 'speech2.wav';

noisefilename = cell(N_noise,3);
%noisefilename{1,1} = 'Babble_noise1.wav';

% length of the recorded microphone signals in seconds
length_recmicsig = 10; %sec

NrOfMics = size(m_pos,1);
NrOfSamples = fs_RIR*length_recmicsig;

speechMatrix = zeros(NrOfSamples,N_speech);
noiseMatrix = zeros(NrOfSamples,N_noise);
mic = zeros(NrOfSamples,NrOfMics);

% read in all speech and noise audiofiles
% resample of speech and noise files and cut off to number of samples
for i =1:N_speech
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:NrOfSamples);
end

for i =1:N_noise
    [noisefilename{i,2}, noisefilename{i,3}] = audioread(noisefilename{i,1});
    noiseTmp = resample(noisefilename{i,2},fs_RIR,noisefilename{i,3});
    noiseMatrix(:,i) = noiseTmp(1:NrOfSamples);
end

% filteroperation
for i=1:NrOfMics
    % filter speech
    for j = 1:N_speech
        mic(:,i) = mic(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
    % add noise
    for j = 1:N_noise
        mic(:,i) = mic(:,i) + fftfilt(RIR_noise(:,i,j),noiseMatrix(:,j));
    end
end

% save files;
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

%% Short Time Fourier Transform
window = 1024;
noverlap = window/2;
nfft = 1024;
NrOfFrames = fix((length(mic(:,1))-noverlap)/(window-noverlap));
NrOfBins = nfft/2 + 1;
S = zeros(NrOfBins,NrOfMics,NrOfFrames);
% F = zeros(NrOfBins,NrOfFrames,NrOfMics);
% T = zeros(NrOfBins,NrOfFrames,NrOfMics);
P = zeros(NrOfBins,NrOfMics,NrOfFrames);

for i = 1:NrOfMics
    x = mic(:,i);
    [S(:,i,:),F,~,P(:,i,:)]= spectrogram(x,window,noverlap,nfft,fs_RIR);
end

% alternative way of calculating power in each frequency bin
% P_total = zeros(NrOfBins,NrOfMics);
% %average out time
% for i = 1:NrOfMics
%     for j =1:NrOfFrames
%      P_total(:,i) = P_total(:,i) + P(:,i,j);
%     end
% end
% %average out mics
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
binCorrMatrix = zeros(NrOfMics, NrOfMics, NrOfBins);
for t=1:NrOfFrames
    for i = 1:NrOfBins
        binCorrMatrix(:,:,i) = binCorrMatrix(:,:,i) + S(i,:,t)'*S(i,:,t);
    end
end

PSD = zeros(NrOfBins,1);
for i = 1:NrOfBins
    PSD(i) = trace(binCorrMatrix(:,:,i));
end

[maxPSD,maxFreqBin] = max(abs(PSD));
maxOmega = 2*pi*F(maxFreqBin);
figure; plot(abs(PSD));
title('Signal power in each frequency bin');xlabel('frequency bins');ylabel('signal power');
hold on; stem(maxFreqBin,maxPSD,'r');

[V,D] = eigs(binCorrMatrix(:,:,maxFreqBin),NrOfMics);
E = V(:,N_speech+1:end);

% dm
dm = m_pos(:,2)-m_pos(1,2);


theta = 0:0.5:180;
g_theta = zeros(NrOfMics, length(theta));
c=340;
for i=1:length(theta)
    g_theta(:,i) = exp( -1i.*maxOmega*dm.*cos(theta(i))./c );
end

% Calculate pseudospectrum for each theta for the maximal freq bin
Pseudospectrum = zeros(length(theta),1);
for i=1:length(theta)
    Pseudospectrum(i) =  1/( (g_theta(:,i))'*E*(E')*g_theta(:,i) );
end

[maxPseudo,maxThetaBin] = max(abs(Pseudospectrum));
DOA_est = theta(maxThetaBin);
figure; plot(0:0.5:180,abs(Pseudospectrum));
title('Pseudospectrum in function of angle theta');xlabel('theta');ylabel('pseudospectrum');
hold on; stem(DOA_est,maxPseudo,'r');

% save estimate of DOA in file
savefile = 'DOA_est.mat';
save(savefile, 'DOA_est');
