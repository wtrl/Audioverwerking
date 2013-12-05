clc;
%% Create the microphone signals
create_micsigs;

Q = N_speech + N_noise;

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


%% delay and sum filter

delay_m = round(-dm*cos((pi/180)*DOA_est(1))*fs_RIR./c);
maxDelay = max(abs(delay_m));
DAS_noise = zeros(nrOfSamples+maxDelay,1);
DAS_speech = zeros(nrOfSamples+maxDelay,1);
if(delay_m(2)>0)
    %noise
    for i = 1:nrOfMics
        DAS_noise = DAS_noise + [zeros(maxDelay-delay_m(i),1); noise(: ,i); zeros(delay_m(i),1)] ;
    end
    %speech
    for i = 1:nrOfMics
        DAS_speech = DAS_speech + [zeros(maxDelay-delay_m(i),1); speech(: ,i); zeros(delay_m(i),1)] ;
    end
else
    %noise
    for i = 1:nrOfMics
        DAS_noise = DAS_noise + [zeros(-delay_m(i),1); noise(: ,i); zeros(maxDelay+delay_m(i),1)] ;
    end
    %speech
    for i = 1:nrOfMics
        DAS_speech = DAS_speech + [zeros(-delay_m(i),1); speech(: ,i); zeros(maxDelay+delay_m(i),1)] ;
    end
end


%speech + noise
DAS_out = DAS_speech + DAS_noise;

%calculate SNR_out

VAD=abs(DAS_speech(:,1))>std(DAS_speech(:,1))*1e-3;
P_speech_DAS = 1./length(VAD==1)*(DAS_speech(VAD==1,1).'*DAS_speech(VAD==1,1));
P_noise_DAS = sum(DAS_noise(:,1).^2)./(nrOfSamples+maxDelay);

SNR_out_DAS = 10*log10(P_speech_DAS./P_noise_DAS)

figure; plot(1:nrOfSamples,mic(:,1),'b',1:nrOfSamples,DAS_out(1:nrOfSamples)*max(mic(:,1))/max(DAS_out),'r');