clear;
clc;
%load HRTFs
load('HRTF.mat');

%resample frequency
fs_RIR = 8000;
N_speech = 2;

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';
speechfilename{2,1} = 'speech2.wav';

% noisefilename = cell(N_noise,3);
% noisefilename{1,1} = 'Babble_noise1.wav'; %'speech2.wav'; 'White_noise1.wav'; 'Babble_noise1.wav'; 

% length of the recorded microphone signals in seconds
length_recmicsig = 10;
% number of samples in the recording length
nrOfSamples = fs_RIR*length_recmicsig;

% Create the microphone signals
speechMatrix = zeros(nrOfSamples,N_speech);
% noiseMatrix = zeros(nrOfSamples,N_noise);

% read in all speech and noise audiofiles
% resampling of speech and noise files and cut off to number of samples
for i = 1:N_speech
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:nrOfSamples);
end

x = speechMatrix(:,1);

binaural_sig1_1 = [x x];
binaural_sig1_2 = [x 0.5*x];
binaural_sig1_3 = [x [0; 0; 0; x(1:end-3)] ];
binaural_sig1_4 = [fftfilt(HRTF(:,1),x) fftfilt(HRTF(:,2),x)];

x = speechMatrix(:,2);

binaural_sig2_1 = [x x];
binaural_sig2_2 = [0.5*x x];
binaural_sig2_3 = [[0; 0; 0; x(1:end-3)] x];
binaural_sig2_4 = [fftfilt(HRTF(:,2),x) fftfilt(HRTF(:,1),x)];

%add them
binaural_sig1 = binaural_sig1_1 + binaural_sig2_1;
binaural_sig2 = binaural_sig1_2 + binaural_sig2_2;
binaural_sig3 = binaural_sig1_3 + binaural_sig2_3;
binaural_sig4 = binaural_sig1_4 + binaural_sig2_4;