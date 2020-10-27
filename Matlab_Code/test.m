% 
% Phase=angle(hilbert(e11_delta));
% 
% Amp=abs(hilbert(e11_gamma)); 
srate=500;
flow1=0.1;
fhigh1=4;
bwlow=25;
bwhigh=40;
[mod2d_norm, p_val_mtx, flow, fhigh] = modulationindex_stats(e11, srate,flow1,fhigh1,bwlow,bwhigh)