%% DAS
%fixed filter = bandpass filter tuned on DOA
DAS_BF;

%% GSC

%filter length of fq
L=1024;
delay = L/2;
mu=0.1;
%Delayed microphone signals
D_mic = D_speech + D_noise;

%Block matrix for case M=5
BM = [1 -1 0 0 0;
      1 0 -1 0 0;
      1 0 0 -1 0;
      1 0 0 0 -1  ];
  
 %Adaptive filter operation
  nrOfDelayedSamples = nrOfSamples + maxDelay;
  noiseRef = zeros(nrOfDelayedSamples+L/2,nrOfMics-1); %In order to be able to listen to the signal out of blocking matrix
  x_k = zeros(L,nrOfMics-1);
  w_k = zeros((nrOfMics-1)*L,1);
  GSC_out = zeros(nrOfDelayedSamples+L/2,1);
  DAS_out = [zeros(L/2,1); DAS_out]; %introduce delay of L/2 taps
  D_mic = [D_mic; zeros(L/2, nrOfMics)];
  VAD = [zeros(L/2,1); VAD];
  for k=1:nrOfDelayedSamples+L/2
      %input to multi channel adaptive filter
      x_k(2:end,:) = x_k(1:end-1,:);
      x_k(1,:) = BM*D_mic(k,:).';
      noiseRef(k,:) = x_k(1,:);
      %filter
      d_pred_k = (w_k.')*x_k(:);
      
      %fixed filter
      d_k = DAS_out(k);
      
      %GSC output
      GSC_out(k) = d_k-d_pred_k;
      
      %update adaptive filtertaps 
      if(VAD(k)==0) %No speech activity
        w_k = w_k + mu./norm(x_k(:)).*x_k(:).*(d_k - d_pred_k); 
      end
  end
  GSC_out = GSC_out(L/2+1:L/2+nrOfDelayedSamples); %remove delay L/2 again
  VAD = VAD(L/2+1:end);

%SNR GSC OUT
P_speech_GSC = 1./length(VAD==1)*(GSC_out(VAD==1).'*GSC_out(VAD==1));
P_noise_GSC = 1./length(VAD==0)*(GSC_out(VAD==0).'*GSC_out(VAD==0));

SNR_out_GSC = 10*log10(P_speech_GSC./P_noise_GSC)

%plots
plot(1:nrOfSamples,GSC_out(plot_offset+1:plot_offset+nrOfSamples),'g',1:nrOfSamples,DAS_speech(plot_offset+1:nrOfSamples+plot_offset)./nrOfMics,'k');