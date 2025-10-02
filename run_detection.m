close all;
%important not to clear the workspace
test_type = 'BH_test'; %either 1) 'BH_test', 'Other_test' (for other bird sounfs),or 'Rand_test'
test_number = 5; %from 1 to 5 for BH and Random sounds, 1 to 6 for Other bird sounds

%%
if strcmp(test_type, 'BH_test')  %get test timestamp if its a BH test signal
    timestamp = bh_timestamps{test_number,1};
end
if strcmp(test_type, 'Rand_test') %to decide the file format
    file = append(test_type, '/', test_type, num2str(test_number), '.mp3');
else
    file = append(test_type, '/', test_type, num2str(test_number), '.aifc');
end
[y, fso] = audioread(file);
fs = 44100; 
y = y(:,1); %only one channel
y = resample(y, fs, fso); %downsample to 44.1 KHz
%k = 1;
%while y(k) < 0.01
%   k  = k + 1;
%end
%y  = y(k:end);
y = spectralsub(y,fs); %spectral subtraction        
y = highpass(y, 6000, fs); %highpass filter
figure
spectrogram(y,[0, floor(length(y)/fs)],fs); %plotting processed test sample spectrogram
%the chosen frame size, overlap etc. is identical to that in
%'make_templates' script
frame_size  = 0.02;
n = 2^nextpow2(fs * frame_size); 
L = floor(0.02* fs);
overlap = floor(L/2); 
step = L - overlap;
tempnum = floor(length(y)/fs);
loop_num = floor((length(y)-fs)/(fs/2));
test_seav = zeros(n/2,loop_num);
test_first = zeros(loop_num,1);
test_second = zeros(loop_num,1);
z = 1;
found = zeros(1,1);
for i = 1:loop_num %looping for each one second frame, to take its seav
        % coefficients
    template = y(1 + (i-1) * fs/2: (i-1) * fs/2 + fs);
    t1_length = length(template);
    frame_num = floor((t1_length - overlap)/(step));
    h_win = hamming(L); %using hamming window for each frame
    frame_coeff = zeros(n/2,frame_num);
    for j = 1:frame_num 
        frame = template((j-1) * step + 1: (j-1) * step + L);
        frame = frame .* h_win;
        fframe = fft(frame,n);
        frame_coeff(:,j) = fframe(1:n/2);
    
    end
    for l = 1:n/2
        test_seav(l,i) = 0;
        for k = 1:frame_num
            test_seav(l,i) = test_seav(l,i) + abs(frame_coeff(l,k));
        end
    end
    for b = 1:n/2 % double the value of coefficients for frequencies above 
        % those of BH sounds, to prevent a possible match
        if  b > 300
            test_seav(b,i) = test_seav(b,i) * 2;
        end
    end
    test_seav(:,i)= test_seav(:,i) ./ max(test_seav(:,i)); %normalizing 



    V = test_seav(:,i) - seav(:,4);
    test_first(i) = sqrt(V' * V); %similarity value of each frame with chirp template
    V = test_seav(:,i) - seav(:,3);
    test_second(i) = sqrt(V' * V); %similarity value with buzz template
    
    if test_first(i) < 2.5 %2.5 is the similarity threshold for the chirp sound
       fprintf('Bee Hummingbird chirp detected at %f sec.\n',(i-1)/2);
       found(z) = (i-1)/2; %store time value where BH sound was detected
       z = z + 1;
    elseif test_second(i) < 1.5
       fprintf('Bee Hummingbird buzz detected at %f sec.\n',(i-1)/2);
       %1.5 for the buzz sound
       found(z) = (i-1)/2;
       z = z + 1;
    end
    
end
if strcmp(test_type, 'BH_test')  %code to check accuracy of BH sound detection,
    %based on given vector pair of BH sound time intervals
    allt = length(timestamp);
    t_count = 0;
    for f = 1:length(timestamp)
        ii = found >= timestamp(f,1) & found <= timestamp(f,2);
        if any(ii)
            t_count = t_count + 1;
        end
    end

    accuracy = t_count / allt;
    fprintf('\nAlgorithm accuracy on BH test sample %d is %d%c \n' , test_number, round(accuracy*100), '%');
end
    
%     figure
%     plot(abs(test_seav(:,77)))
%     figure
%     plot(abs(test_seav(:,76)))
    %plotting the two different test template seav vectors
    figure
    plot(seav(:,4)); 
    figure
    plot(seav(:,3));
