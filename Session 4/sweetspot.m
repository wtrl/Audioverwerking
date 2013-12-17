clear;
clc;
%load HRTFs
load('HRTF.mat');
load('Computed_RIRs.mat');
nrOfLoudspeakers = size(s_pos,1);
Lh = 1500;
M = nrOfLoudspeakers;
Lg = ceil(2*(Lh-1)/(M-2));
%pre-defined delay
Delta = ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);

%HRTF
% xL = 1;
% xR = 1;
nrOfBinauralSigs = 4;
xL = [1; zeros(1000,1)];
binaural_sig1_1 = [xL xL];
binaural_sig1_2 = [xL 0.5*xL];
binaural_sig1_3 = [xL [0; 0; 0; xL(1:end-3)] ];
binaural_sig1_4 = [fftfilt(HRTF(:,1),xL) fftfilt(HRTF(:,2),xL)];

x_binaural_sig = cell(nrOfBinauralSigs,1);
x_binaural_sig{1} = binaural_sig1_1 ;
x_binaural_sig{2} = binaural_sig1_2 ;
x_binaural_sig{3} = binaural_sig1_3 ;
x_binaural_sig{4} = binaural_sig1_4 ;
%System-Of-Equations
%H
RIR_sources_truncated = RIR_sources(1:Lh,:,:);
H_L = zeros(Lh+Lg-1,M*Lg);
H_R = zeros(Lh+Lg-1,M*Lg);
for i=1:M
    H_L(:, 1+(i-1)*Lg : Lg*i) = toeplitz([RIR_sources_truncated(:,1,i); zeros(Lg-1,1)],[RIR_sources_truncated(1,1,i); zeros(Lg-1,1)]);
    H_R(:, 1+(i-1)*Lg : Lg*i) = toeplitz([RIR_sources_truncated(:,4,i); zeros(Lg-1,1)],[RIR_sources_truncated(1,2,i); zeros(Lg-1,1)]);
end
H = [H_L; H_R];
H_L_real = zeros(Lh+Lg-1,M*Lg);
H_R_real = zeros(Lh+Lg-1,M*Lg);
for i=1:M
    H_L_real(:, 1+(i-1)*Lg : Lg*i) = toeplitz([RIR_sources_truncated(:,2,i); zeros(Lg-1,1)],[RIR_sources_truncated(1,1,i); zeros(Lg-1,1)]);
    H_R_real(:, 1+(i-1)*Lg : Lg*i) = toeplitz([RIR_sources_truncated(:,5,i); zeros(Lg-1,1)],[RIR_sources_truncated(1,2,i); zeros(Lg-1,1)]);
end
H_real = [H_L_real; H_R_real];
%remove zero rows
ix = any(H,2);
H_small = H(ix,:);
g = zeros(M*Lg,nrOfBinauralSigs);
for j = 1:nrOfBinauralSigs
    
    %x
    xL = x_binaural_sig{j}(:,1);
    xR = x_binaural_sig{j}(:,2);
    x = [zeros(Delta,1); xL; zeros(Lg+Lh-1-Delta-length(xL),1); zeros(Delta,1); xR; zeros(Lg+Lh-1-Delta-length(xR),1)];
    x = x(ix);
    
    
    g(:,j) = H_small\x;
    
    figure; plot(1:size(x,1),H_small*g(:,j),'b',1:size(x,1), x,'r');
    synth_error = norm(H_small*g(:,j)-x)
    
end

%resample frequency
fs_RIR = 8000;
N_speech = 1;

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';
% speechfilename{2,1} = 'speech2.wav';

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

speech = speechMatrix(:,1);

Left=zeros(length(speech),nrOfBinauralSigs);
Right= zeros(length(speech),nrOfBinauralSigs);

for j=1:nrOfBinauralSigs
    for i=1:nrOfLoudspeakers
        Left(:,j) = Left(:,j) + fftfilt(H_L_real(:, 1+(i-1)*Lg : Lg*i)*g( 1+(i-1)*Lg : Lg*i,j), speech);
        Right(:,j) = Right(:,j) + fftfilt(H_R_real(:, 1+(i-1)*Lg : Lg*i)*g( 1+(i-1)*Lg : Lg*i,j), speech);
    end
end

binaural_sig1 = [Left(:,1) Right(:,1)];
binaural_sig2 = [Left(:,2) Right(:,2)];
binaural_sig3 = [Left(:,3) Right(:,3)];
binaural_sig4 = [Left(:,4) Right(:,4)];
