clear;

% user defined file names of the WAV files to use for sources and noise
load('Computed_RIRs.mat');

%Number of files
N_speech = size(s_pos,1);
N_noise = size(v_pos,1);

speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';
%speechfilename{2} = 'speech2.wav';

noisefilename = cell(N_noise,3);
%noisefilename{1,1} = 'Babble_noise1.wav';

% length of the recorded microphone signals in seconds
length_recmicsig = 5;

NrOfMics = size(m_pos,1);
NrOfSamples = fs_RIR*length_recmicsig;

speechMatrix = zeros(NrOfSamples,N_speech);
noiseMatrix = zeros(NrOfSamples,N_noise);
mic = zeros(NrOfSamples,NrOfMics);

% read in all speech and noise audiofiles
for i =1:N_speech
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
end

for i =1:N_noise
    [noisefilename{i,2}, noisefilename{i,3}] = audioread(noisefilename{i,1});
end

% resample of speech and noise files and cut off to number of samples
for i =1:N_speech
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:NrOfSamples);
end

for i =1:N_noise
    noiseTmp = resample(noisefilename{i,2},fs_RIR,noisefilename{i,3});
    noiseMatrix(:,i) = noiseTmp(1:NrOfSamples);
end

%filteroperation
for i=1:NrOfMics
    %filter speech
    for j = 1:N_speech
        mic(:,i) = mic(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
    %add noise
    for j = 1:N_noise
        mic(:,i) = mic(:,i) + fftfilt(RIR_noise(:,i,j),noiseMatrix(:,j));
    end
end


%save files;
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

figure; plot(mic(:,1),'r');
hold on
plot(mic(:,2),'b');
hold off