%all_callbacks Contains all the callbacks of mySA_GUI
%
% Syntax:  all_callbacks()
%
% Inputs:
%
% Outputs:
%
%
% Example: 
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
% MATLAB version: R2012b
%
% See also: 
%
% Authors: Giuliano Bernardi
% KU Leuven, Department of Electrical Engineering (ESAT/SCD)
% email: giuliano.bernardi@esat.kuleuven.be
% Website: http://homes.esat.kuleuven.be/~gbernard/index.html
% Created: 22-October-2013; Last revision: 25-October-2013

%------------- BEGIN CODE --------------

set(S.eh_save,'callback',{@mh_save_session,S});
set(S.eh_load,'callback',{@mh_load_session,S});
set(S.eh_help,'callback',{@mh_help_about,S});
set(S.eh_about,'callback',{@mh_help_about,S});
set(S.ax,'ButtonDownFcn',{@ax_bdfcn,S});
set(S.pb_audio,'callback',{@pb_call,S});  
set(S.pb_audio_rm,'callback',{@pb_del_single_cmp,S});  
set(S.pb_noise,'callback',{@pb_call,S});
set(S.pb_noise_rm,'callback',{@pb_del_single_cmp,S});  
set(S.pb_mic,'callback',{@pb_call,S});
set(S.pb_mic_rm,'callback',{@pb_del_single_cmp,S});  
set(S.pb_reset,'callback',{@pb_reset,S});
set(S.pb_create_RIRs,'callback',{@pb_create_RIRs,S});
set(S.ed_nmics,'callback',{@ed_kpfcn,S});
set(S.ed_dmics,'callback',{@ed_kpfcn,S});
set(S.ed_rdim,'callback',{@ed_par_kpfcn,S});
set(S.ed_reverb,'callback',{@ed_par_kpfcn,S});
set(S.ed_lRIR,'callback',{@ed_par_kpfcn,S});
set(S.pb_DOA,'callback',{@pb_draw_DOA,S});
set(S.pb_reset_DOA,'callback',{@pb_reset_DOAs,S});
set(S.pb_reset_cmp,'callback',{@pb_reset_cmp,S});
set(S.pp_fs,'callback',{@pp_get_fs,S});

%------------- END CODE --------------