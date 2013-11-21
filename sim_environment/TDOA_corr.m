clear;

load('Computed_RIRs.mat');

NrOfMics = size(m_pos,1);

indicesMax = zeros(NrOfMics,1);

for i=1:NrOfMics
    [~,indicesMax(i)] = max(RIR_sources(:,i));
end

ground_truth_TDOA = indicesMax(1)-indicesMax(2)


%Number of files
N_speech = size(s_pos,1);
N_noise = size(v_pos,1);

% user defined file names of the WAV files to use for sources and noise
speechfilename = cell(N_speech,3);
speechfilename{1,1} = 'speech1.wav';

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
%     speechMatrix = wgn(NrOfSamples,1,1);
    [speechfilename{i,2}, speechfilename{i,3}] = audioread(speechfilename{i,1});
    speechTmp = resample(speechfilename{i,2},fs_RIR,speechfilename{i,3});
    speechMatrix(:,i) = speechTmp(1:NrOfSamples);
end

% for i =1:N_noise
%     [noisefilename{i,2}, noisefilename{i,3}] = audioread(noisefilename{i,1});
% end


% for i =1:N_noise
%     noiseTmp = resample(noisefilename{i,2},fs_RIR,noisefilename{i,3});
%     noiseMatrix(:,i) = noiseTmp(1:NrOfSamples);
% end

%filteroperation
for i=1:NrOfMics
    %filter speech
    for j = 1:N_speech
        mic(:,i) = mic(:,i) + fftfilt(RIR_sources(:,i,j),speechMatrix(:,j));
    end
    %add noise
%     for j = 1:N_noise
%         mic(:,i) = mic(:,i) + fftfilt(RIR_noise(:,i,j),noiseMatrix(:,j));
%     end
end


%save files;
savefile = 'mic.mat';
save(savefile, 'mic', 'fs_RIR');

figure; plot(mic(:,1),'r'); title('The 2 microphone signals');
hold on
plot(mic(:,2),'b');
hold off

% cross correlation method
mic1_crossc = mic(:,1);%mic(7300:26000,1);
figure; plot(mic1_crossc); title('Representative signal segment');
c = xcorr(mic1_crossc,mic(:,2));
% c = c(NrOfSamples:end);
[max_c,estimated_TDOA] = max(c);
estimated_TDOA = estimated_TDOA-NrOfSamples
figure; plot(c,'b'); title('Time domain cross-correlation function with ground truth estimate');
hold on; stem(NrOfSamples+ground_truth_TDOA,max_c,'r');

diff_TDOA = abs(estimated_TDOA-ground_truth_TDOA)