TDOA_corr;

% speed of sound in air in m/s
c = 340;

% inter-microphone distance in m
d = norm(m_pos(1,:)-m_pos(2,:));

% DOA estimate in degrees between 0-180 degrees
arg = estimated_TDOA*c/(d*fs_RIR);
if(arg > 1)
    DOA_est = 0;
elseif(arg < -1)
    DOA_est = 180;
else
    DOA_est = acos(estimated_TDOA*c/(d*fs_RIR))*180/pi;
end
        
%save estimate of DOA in file
savefile = 'DOA_est.mat';
save(savefile, 'DOA_est');