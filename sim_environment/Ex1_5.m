clear;

load('Computed_RIRs.mat');

NrOfMics = size(m_pos,1);

%Number of files
N_speech = size(s_pos,1);
N_noise = size(v_pos,1);

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';
speechfilename{2,1} = 'speech2.wav';

% length of the recorded microphone signals in seconds
length_recmicsig = 5;

NrOfSamples = fs_RIR*length_recmicsig;

speechMatrix = zeros(NrOfSamples,N_speech);

mic = zeros(NrOfSamples,NrOfMics);

% read in all speech and noise audiofiles
for i = 1:N_speech
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:NrOfSamples);
end

% speechMatrix = wgn(NrOfSamples,N_speech,1);

%filteroperation
for i = 1:NrOfMics
    %filter speech
    for j = 1:N_speech
        mic(:,i) = mic(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
end

%save files;
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

figure; plot(mic(:,1),'r'); title('The 2 microphone signals');
hold on
plot(mic(:,2),'b');
hold off

% cross correlation method
c = xcorr(mic(:,1),mic(:,2));
[max_c,estimated_TDOA] = max(c);
estimated_TDOA = estimated_TDOA-NrOfSamples;
figure; plot(c,'b'); title('Time domain cross-correlation function');