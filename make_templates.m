close all;
%important not to clear the workspace
fs = 44100; %sampling rate, no need to downsample to 44.1 KHz, 
            %already done in run_identifier script
frame_size = 0.02; %frame size in sec
overlap = 0.01; %overlap is 
to  = 0;
tf = 0;
n = 2^nextpow2(fs * frame_size);
seav = zeros(n/2,4); %matrix that will contain seav vector of all 4 templates 
for g = 3:6 %only analyzing 4 last signals from BH_database, which contain instances
    %4 templates, 2 distinct sound patterns, 2 templates for each sound
    %pattern
    %of perfect BH sound patters
    if g == 3 || g == 6 %for each signal, take second where distinct pattern appears
        to = 2; 
        tf = 3;
    elseif g == 4
        to = 3;
        tf = 4;
    else 
        to = 1.4;
        tf = 2.4;
    end    
    y = sample{g,1}; %sample 4 from which template is taken  
    t1 = y(floor(to*fs):floor(tf*fs)); %first template vector
    t1_length = length(t1);
    L = floor(fs*0.02); %frame length in samples
    overlap = floor(L/2); %half the frame size
    step = L - overlap; %step for every frame
    frame_num = floor((t1_length - overlap)/(step)); %number of frames in template sec.
    h_win = hamming(L); %using hamming window for each frame
    frame_coeff = zeros(n/2,frame_num);
    for i = 1:frame_num %loop to get n/2 fft coefficents for all frames of template
        frame = t1((i-1) * step + 1: (i-1) * step + L);
        frame = frame .* h_win;
        fframe = fft(frame,n);
        frame_coeff(:,i) = fframe(1:n/2); %stored in this (n/2)-by-frame_num matrix
    
    end
    for i = 1:n/2
        seav(i,g-2) = 0; %matrix where seav vectors for each 1-sec template sound are stored
        for k = 1:frame_num
            seav(i,g-2) = seav(i,g-2) + abs(frame_coeff(i,k));
        end
    end
    seav(:,g-2)= seav(:,g-2) ./ max(seav(:,g-2)); %normalizing the seav vectors
    figure
    plot(abs(seav(:,g-2)));
end
V = seav(:,1) - seav(:,4);
d_first = sqrt(V' * V);
V = seav(:,2) - seav(:,3);
d_second = sqrt(V' * V);
%d_first, d_second are values to give idea on where the similarity threshold
%should be placed, by comparing the pair of seav vectors corresponding
%to the same 1-sec sound pattern 