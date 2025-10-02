close all;
clear all;
samplenum = 6; %number of samples in BH database path
prefix = 'mellisuga/BH_Sample'; %name of the BH database path
suffix1 = '.aifc'; %file format
sample = cell(samplenum,1); %
cfs = 44100; %the standard sampling frequency used to downsample the signals
for i = 1:samplenum
    s_name = append(prefix,int2str(i));
    s_name = append(s_name, suffix1);
    
    [y,fs] = audioread(s_name);
    y = y(:,1);
    y = resample(y, cfs, fs); %downsampling each signal to 44.1 KHz
    k = 1;
    while y(k) < 0.01 %loop to remove total silence from beginning of signal
        k  = k + 1;
    end
    y  = y(k:end);
    y = spectralsub(y,cfs); %signal subtraction noise cancelling on signal
        
    y = highpass(y, 6000, cfs); %highpass filter with 6 KHz cut-off
    sample{i,1} = y;
    figure;
    spectrogram(y,[0, floor(length(y)/cfs)],cfs); %Plotting signal Spectrogram
    titl = append(prefix, int2str(i));
    title(titl);
    figure ;
    plot(y); %plotting signal
    title(titl);
    saveas(gcf,append(titl,'.png')); %saving signal image

end


