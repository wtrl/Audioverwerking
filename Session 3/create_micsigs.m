clear;

% user defined file names of the WAV files to use for sources and noise
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

noisefilename = cell(N_noise,3);
noisefilename{1,1} = 'speech2.wav'; %'speech2.wav'; 'White_noise1.wav'; 'Babble_noise1.wav'; 

% length of the recorded microphone signals in seconds
length_recmicsig = 10;
% number of samples in the recording length
nrOfSamples = fs_RIR*length_recmicsig;

%% Create the microphone signals
speechMatrix = zeros(nrOfSamples,N_speech);
noiseMatrix = zeros(nrOfSamples,N_noise);

% read in all speech and noise audiofiles
% resampling of speech and noise files and cut off to number of samples
for i = 1:N_speech
    [speechfilename{i,2}, speechfilename{i,3}] = wavread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:nrOfSamples);
end

for i = 1:N_noise
    [noisefilename{i,2}, noisefilename{i,3}] = wavread(noisefilename{i,1});
    noiseTmp = resample(noisefilename{i,2},fs_RIR,noisefilename{i,3});
    noiseMatrix(:,i) = noiseTmp(1:nrOfSamples);
end
%noiseMatrix = wgn(nrOfSamples,N_noise,1,'dBm'); % use wgn instead of noise signal

speech = zeros(nrOfSamples,nrOfMics);
noise = zeros(nrOfSamples,nrOfMics);
% filter operation
for i = 1:nrOfMics
    % filter speech
    for j = 1:N_speech
        speech(:,i) = speech(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
    
    for j = 1:N_noise
        noise(:,i) = noise(:,i) + fftfilt(RIR_noise(:,i,j),noiseMatrix(:,j));
    end
    
end
power = 10*log10(sum(speech(:,1).^2)./nrOfSamples);
% add white gaussian noise
noise = wgn(nrOfSamples,nrOfMics,power-10,'dbW')+noise;
mic = speech +noise;

% calculate noise covariance matrix
% noiseCovarianceMatrix = noiseOnMic'*noiseOnMic;

VAD=abs(speech(:,1))>std(speech(:,1))*1e-3;

P_speech = 1./length(VAD==1)*(speech(VAD==1,1).'*speech(VAD==1,1));
P_noise = sum(noise(:,1).^2)./nrOfSamples;

SNR_in = 10*log10(P_speech./P_noise)
% save microphone signals and sampling rate in file
savefile = 'mic.mat';
save(savefile, 'mic','speech', 'noise', 'SNR_in', 'fs_RIR');
