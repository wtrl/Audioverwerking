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
  
  %input to multi channel adaptive filter
  nrOfSamples = nrOfSamples + maxDelay;
  noiseRef = zeros(nrOfSamples,nrOfMics-1);
  for k=1:nrOfSamples
      noiseRef(k,:) = BM*D_mic(k,:).';
  end
  noiseRef = [zeros(L-1,nrOfMics-1); noiseRef;zeros(L/2,nrOfMics-1);];

%   w_k = zeros((nrOfMics-1)*L,nrOfSamples+1);
  w_k = zeros((nrOfMics-1)*L,1);
  GSC_out = zeros(nrOfSamples+L/2,1);
  DAS_out = [zeros(L/2,1); DAS_out]; %introduce delay of L/2 taps
  for k=1:nrOfSamples+L/2
      %adaptive filter
      x_k = noiseRef(k:k+L-1,:);
      d_pred_k = (w_k.')*x_k(:);
      
      %fixed filter
      d_k = DAS_out(k);
      
      %GSC output
      GSC_out(k) = d_k-d_pred_k;
      
      %update adaptive filtertaps
%       w_k(:,k+1) = w_k(:,k) + mu./norm(x(:)).*x(:).*(d_k - d_pred_k); 
      w_k = w_k + mu./norm(x_k(:)).*x_k(:).*(d_k - d_pred_k); 

  end
  GSC_out = GSC_out(L/2+1:L/2+nrOfSamples);

P_speech_GSC = 1./length(VAD==1)*(GSC_out(VAD==1).'*GSC_out(VAD==1));
P_noise_GSC = 1./length(VAD==0)*(GSC_out(VAD==0).'*GSC_out(VAD==0));

SNR_out_GSC = 10*log10(P_speech_GSC./P_noise_GSC)

%plots
if(delay_m(2)>0)
    plot(1:nrOfSamples-maxDelay,GSC_out(maxDelay+1:nrOfSamples),'g',1:nrOfSamples,DAS_speech(maxDelay+1:nrOfSamples)./nrOfMics,'k');
else
    plot(1:nrOfSamples-maxDelay,GSC_out(1:nrOfSamples-maxDelay),'g',1:nrOfSamples,DAS_speech(1:nrOfSamples-maxDelay)./nrOfMics,'k');
end