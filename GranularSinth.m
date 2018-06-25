clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%synthesis parameters and options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs_sys = 48000; %System sampling rate Hz
filePath = 'Voice.wav'; % sound file
grainSelMode = 2; % 1 = fixed position + spray; 2 = sequential + spray; 3 = random
duration_sec = 4; %length of output signal in seconds (any value above 0s)
grainSize_msec = 75; %size of grain in milli seconds (any value above 0s and below the lenght of the sound file)
grainOverlap_perc = 0.0; % percentage of overlap between consecutive grains (any value between 0 and 1 (not included))  
kaiserBeta = 0; %defines the window shape, 0 -> retangular, ?? --> Hanning, ?? --> Hamming
fftSize_smp = 512; % size fo FFT for spectrum computation, must be power of 2
phaseRateScale = 5; % phase adjust rate scaling factor;

filePos_perc = 0.3; % file position - percentage between 0 and 1 - to extract grains (mode 1 and 2 only)
spray_perc = 0.7; % random deviation in percentage between 0 and 1 from the fixed position to extract grains (mode 1 and 2 only)
spray_mode =3; % 1 backward and forward, 2 backward only, 3 forward only (mode 1 and 2 only)
sequenceRate = 0.3; % grain selection rate (positive or negative) with respect grain size (mode 2 only)


%%parameter check, conversion and sound loading

[snd,Fs_sig]=audioread(filePath); %read audio
%sound(snd,Fs_sig) 

%check number of channels -> change to 1
if size(snd,2) > 1
    snd = snd (:,1); %if it is multiple then use the first channel
end

%change sampling rate to 48k
if Fs_sig ~= Fs_sys
    snd = resample(snd, Fs_sys, Fs_sig);
end 

sndLength_smp = length(snd); %length of input signal in samples

if duration_sec < 0
    error('negative duration');
end

if (grainSize_msec < 0) || (grainSize_msec/1000 > (sndLength_smp/Fs_sys))
    error('grain size out of range');
end

if (sum(grainSelMode == [1,2,3])==0)
    error('not supported grain selection mode');
end

if (filePos_perc < 0) || (filePos_perc > 1)
    error('file position out of range');
end

if (spray_perc < 0) || (spray_perc > 1)
    error('spray percentage out of range');
end

duration_smp = round(duration_sec * Fs_sys); %length of output signal in samples
grainSize_smp = round((grainSize_msec / 1000) * Fs_sys); %size of grain in samples 
grainStep_smp = round(grainSize_smp - (((grainOverlap_perc * grainSize_msec) / 1000) * Fs_sys)); %grain step size in samples
filePos_smp = round(filePos_perc * (sndLength_smp - grainSize_smp - 1)); %file position in samples for mode 1

if ~((grainStep_smp >= 1) || (grainStep_smp <= grainSize_smp))
    error('invalid step size')
end

if fftSize_smp < grainSize_smp
    fftSize_smp = 2^ceil(log2(grainSize_smp));
    warning('fft size too small, increasing it to match grain size');
end

if (rem(fftSize_smp,2)~=0)
    fftSize_smp = fftSize_smp + 1;
end

if grainSelMode == 1 % fixed + spray
    
    if spray_mode == 1 % 1 backward and forward, 2 backward only, 3 forward only
        if (filePos_smp > ((sndLength_smp - grainSize_smp - 1) - filePos_smp))
            spray_smp = 2 * round(((sndLength_smp - grainSize_smp - 1) - filePos_smp) * spray_perc);
            sprayOffset_smp = spray_smp / 2;
        else
            spray_smp = 2 * round(filePos_smp * spray_perc);
            sprayOffset_smp = spray_smp / 2;
        end
    elseif spray_mode == 2
        spray_smp = - round(filePos_smp * spray_perc);
        sprayOffset_smp = 0;
    elseif spray_mode == 3
        spray_smp = round(((sndLength_smp - grainSize_smp - 1) - filePos_smp) * spray_perc);
        sprayOffset_smp = 0;
    else
        error('wrong spray mode');
    end
    
elseif grainSelMode == 2 % sequential + spray
    
    %grainSize_smp = abs(round(grainSize_smp * (sequenceRate))); % stretching the grain to match the rate
    
    if spray_mode == 1 % 1 backward and forward, 2 backward only, 3 forward only
        spray_smp = 2 * round(abs(round(grainSize_smp * (sequenceRate))) * spray_perc);
        sprayOffset_smp = spray_smp / 2;
    elseif spray_mode == 2 || spray_mode == 3 
        spray_smp = (abs(round(grainSize_smp * (sequenceRate))) * spray_perc * sign(sequenceRate));
        sprayOffset_smp = 0;
    else
        error('wrong spray mode');
    end
    
end


% generating window function
window = kaiser(grainSize_smp, kaiserBeta); %grain window
% figure(1);
% plot(window);
% title('Window function');



% synthesis loop
niter = ceil((duration_smp - grainSize_smp) / grainStep_smp); % number of iterarions
out_t = zeros(duration_smp,1); %initializing output array
out_s = zeros(duration_smp,1); %initializing output array
phase = zeros(fftSize_smp/2,1);%inizializing phase temporary array

for i=0:(niter - 1)
    
    if grainSelMode == 1
        grain = window .* getgrain_fix(snd, grainSize_smp, filePos_smp, spray_smp, sprayOffset_smp); % get a grain and apply window
    elseif grainSelMode == 2
        grain = window .* getgrain_seq(snd, grainSize_smp, filePos_smp, spray_smp, sprayOffset_smp, sequenceRate, sndLength_smp, i);
    else
        grain = window .* getgrain_rnd(snd, grainSize_smp); % get a grain and apply window
    end
    %time domain granula synthesis
    out_idx = (i * grainStep_smp ) + 1; % index of current grain in output array
    out_t(out_idx:(out_idx + grainSize_smp - 1)) = out_t(out_idx:(out_idx + grainSize_smp - 1)) + grain; %add and overlap current grain to output
    
    %spectral domain granular synthesis
    magspect = abs(fftR2wPad(grain,fftSize_smp));
    magspect = magspect(1:fftSize_smp/2);
    [recon, phase] = SpectrumInversion (magspect, phase, fftSize_smp, grainSize_smp, grainStep_smp, phaseRateScale);
    recon = real(ifft(recon,fftSize_smp));
    recon = recon(1:grainSize_smp).*window;
    out_s(out_idx:(out_idx + grainSize_smp - 1)) = out_s(out_idx:(out_idx + grainSize_smp - 1)) + recon; %add and overlap current grain to output
    
end

sound(out_s,Fs_sys);

out_t = out_t/(max(abs(out_t))); %normalizing amplitude to unitary range
out_s = out_s/(max(abs(out_s))); %normalizing amplitude to unitary range
figure(2);
plot(out_t);
title('time domain granular synthesis');
figure(3);
plot(out_s);
title('frequency domain granular synthesis');

N = 1024;
[s1,f1,t1] = spectrogram(out_t,hamming(N),ceil(0.75*N),N,Fs_sys,'yaxis');
s1 = abs(s1);
sizet = (size(s1,2)-1);
for i = 1:sizet
     diff1 = abs(s1(:,i) - s1(:,i+1)); %positive
     diffMean1(i) = mean(diff1);
end

t1 = t1(1:sizet);
     
figure(4);
plot(t1,diffMean1)
xlabel('time (s)')
ylabel('mean of difference')
title('time domain approach')

[s2,f2,t2] = spectrogram(out_s,hamming(N),ceil(0.75*N),N,Fs_sys,'yaxis');
s2 = abs(s2);
for i = 1:sizet
     diff2 = abs(s2(:,i) - s2(:,i+1)); %positive
     diffMean2(i) = mean(diff2);
end

t2 = t2(1:sizet);
     
figure(5);
plot(t2,diffMean2)
xlabel('time (s)')
ylabel('mean of difference')
title('frequency domain approach')

[pks,locs] = findpeaks(diffMean1,t1,'MinPeakDistance',round((floor(grainSize_smp / (N - ceil(0.75*N)))))/Fs_sys);
[pks2,locs2] = findpeaks(diffMean2,t2,'MinPeakDistance',round((floor(grainSize_smp / (N - ceil(0.75*N)))))/Fs_sys);

sizep1 = size(locs(1,:));
addPeaks1 = 0;
for i = 1:sizep1(1,2)
    addPeaks1 = pks(i)+ addPeaks1;
end 

peaksAvg1 = addPeaks1/sizep1(1,2);

sizep2 = size(locs2(1,:));
addPeaks2 = 0;
for i = 1:sizep2(1,2)
    addPeaks2 = pks2(i)+ addPeaks2;
end 

peaksAvg2 = addPeaks2/sizep2(1,2);

% max(pks)/2
% count = 0
% for i = 1:sizep1(1,2)
%     if pks[i]>(max(pks)/2)
%         count = count +1;
%     end 
% end 

