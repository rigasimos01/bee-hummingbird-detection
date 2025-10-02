function output=spectralsub(signal,fs)

%returns a "clean" version of the input audio signal
% inputs: signal- the data vector of the audio signal to be de-noised
%         fs - the sampling ffrequency of the audio signal


n_size= 0.5; %size in seconds of initital noise segment in audio 

L = floor(0.02 * fs); %frame size
h_win = hamming(L); %hamming window
n = L; %n-point fft resolution
step = floor(0.4 * L); 
overlap = L - step;
nframe = floor((length(signal) - overlap)/step);

ffts = zeros(n, nframe);
mags = zeros(floor(n/2) + 1, nframe);
phases = zeros(floor(n/2) + 1, nframe);
for i = 1:nframe %loop to create magnitutde and phase vectors for all frames, store them to a matrix
    frame = signal((i-1) *step + 1: (i-1) * step + L);
    frame = frame .* h_win;
    fframe = fft(frame,n);
    ffts(:,i) = fframe;
    mags(:,i) = abs(fframe(1:floor(n/2)+1));
    phases(:,i) = angle(fframe(1:floor(n/2)+1));
end

noise_segnum=floor((n_size*fs- L)/(step) +1);%counting frames of initial noise


noise=mean(mags(:,1:noise_segnum)')'; %getting a mean frame for all noise frames
nresmax=zeros(length(noise));% noise residual maximum (taken by SSBol implementation)


mags_avg= mags; 
for i=2:nframe - 1
    mags_avg(:,i)=(mags(:,i-1)+mags(:,i)+mags(:,i+1))/3; %making an average spectrum for each frame, based on previous and next frame
end

for i=1:nframe
    sub=mags_avg(:,i)-noise; % performing ss
    if i>1 && i<nframe %reducing the residual noise after sutraction           
        for j=1:length(sub)
            if sub(j)<nresmax(j)
                    sub(j)=min([sub(j) mags_avg(j,i-1)-noise(j) mags_avg(j,i+1)-noise(j)]);
            end
        end
    end
    clean_mags(:,i)= sub ;
end

recon=clean_mags.*exp(1i*phases);

if mod(L,2) %check if odd
    recon=[recon;flipud(conj(recon(2:end,:)))];
else
    recon=[recon;flipud(conj(recon(2:end-1,:)))];
end

output=zeros((nframe-1)* step+ L,1);
for i=1:nframe %reconstruct clean audio signal 
    f_recon=recon(:,i);
    output((i-1)*step+1:(i-1)*step+L) = output((i-1)*step+1:(i-1)*step+L)+real(ifft(f_recon,L));
end

