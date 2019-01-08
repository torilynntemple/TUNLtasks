function [ms] = msExtractBinary(ms);
%MSEXTRACTBINARY Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
dt = median(diff(ms.time))/1000; % Conversion from ms to s
Fs = 1/dt;
z_threshold = 2;

[bFilt,aFilt] = butter(2,  2/(30/2), 'low');
ms.numNeurons = length(ms.RawTraces(1,:)); 

for trace_i = 1:ms.numNeurons;
raw_trace = ms.RawTraces(:,trace_i);
filt_trace = zscore(filtfilt(bFilt,aFilt,raw_trace));
d1_trace = diff(filt_trace);
d1_trace(end+1) = 0;
d2_trace = diff(d1_trace);
d2_trace(end+1) = 0;

binary_trace = filt_trace*0;
binary_trace(filt_trace>z_threshold & d1_trace>0) = 1;
    
ms.Binary(:,trace_i) = binary_trace;
end



end

