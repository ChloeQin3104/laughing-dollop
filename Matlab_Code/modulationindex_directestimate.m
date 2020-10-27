%   adapted from the code written for Kramer et al. (2008) 
%    by Tolga Ozkurt (2011)
%   to implement direct PAC estimate (Ozkurt and Schnitzler, 2011)
%   and raw moduation index (Canolty et el., 2006)     
%
%   INPUTS:
%    d      = The unfiltered data.
%    srate  = The sampling rate (e.g., 1000 Hz)
%
%   OUTPUS:
%    
%    mod2d  = The two-dim modulation index.
%    flow   = The frequency axis for the lower phase frequencies.
%    fhigh  = The frequency axis for the higher amplitude frequencies.
%


 
function [mod2dDirectEst, mod2dRawMI, flow, fhigh] = modulationindex_directestimate(d, srate,flow1,fhigh1,bwlow,bwhigh)
 
  %The phase frequencies.
  flow2 = flow1+bwlow;                               
  
  %The amp frequencies.
  fhigh2 = fhigh1+bwhigh;                           
  
  mod2dRawMI = zeros(length(flow1), length(fhigh1));
  mod2dDirectEst = zeros(length(flow1), length(fhigh1));
  
  for i=1:length(flow1)
      theta=eegfilt(d,srate,flow1(i),flow2(i));  %Compute the low freq signal.
      theta=theta(srate:length(d)-srate-1);           %Drop the first and last second.
      phase = angle(hilbert(theta));                  %Compute the low freq phase.
      ['Loops remaining = ', num2str(length(flow1)-i)]
      for j=1:length(fhigh1)
        gamma=eegfilt(d,srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
        gamma=gamma(srate:length(d)-srate-1);         %Drop the first and last second.
        amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.

        %Compute the modulation index.
        [m_norm1 m_norm2] = modulation_index_nostats(amp, phase);
        mod2dRawMI(i,j) = m_norm1;
        mod2dDirectEst(i,j) = m_norm2;
      end
  end
 
  %Plot the two-dimensional PAC portraits
  
  
  flow = (flow1 + flow2) / 2;
  fhigh = (fhigh1 + fhigh2) / 2;
  
  figure
  imagesc(flow, fhigh, mod2dRawMI');  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  xlabel('Phase Frequency [Hz]');  ylabel('Amplitude Frequency [Hz]');
  title('Canolty Modulation index')

  figure
  imagesc(flow, fhigh, mod2dDirectEst');  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  xlabel('Phase Frequency [Hz]');  ylabel('Amplitude Frequency [Hz]');
  title('Direct estimator')
  
  % end of the function.
 
  function [m_raw1 m_raw2] = modulation_index_nostats(a, p)
 
  N = length(a);
  z = a.*exp(1i*p);
  m_raw1 = abs(mean(z));   %Compute the mean length 
  m_raw2 = (1./sqrt(N)) * abs(mean(z)) / sqrt(mean(a.*a)); % compute the direct estimate
