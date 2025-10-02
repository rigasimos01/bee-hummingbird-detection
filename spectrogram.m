function imag = spectrogram(signal, time, fs)
% inputs: signal, the vector audio signal, time, 
% a vector with two entries, time of start and time of ending 
% in seconds, and fs, the sampling frequency
% plots narrowband spectorgram of signal
narrowband = 50/8000 * fs;      % window length for narrowband spectrogram
wL = round(1.185 * fs / narrowband);  
h_win = hamming(wL);

db_range = 80;  % dB level range
specslice = 500; %spectrogram slices
n = 2048; %determining frequency resolution 
t = linspace(time(1), time(2), specslice); %define time period of graph

% Get spectrum in dB
t_length = length(t);

freq_db = zeros(n / 2 + 1, t_length);
is = floor(t * fs - (wL-1) / 2) + 1;
%constructing db value matrix, one vector for each slice
for i =1:t_length
    slice = getSlice(signal, h_win, is(i)); 
    fslice = fft(slice, n);
    mags = abs(fslice(1:n / 2 + 1));
    freq_db(:,i) = (20 * log10 (mags));   
end

% Frequency values
freq = ((0:n/2) / n) * fs;

% frequency range
Fmin = 0;
Fmax = fs/2;
ss = find(freq >= Fmin & freq <= Fmax);
freq = freq(ss);
freq_db = freq_db(ss,:);

%gain calibration
Apeak = abs (sum(h_win));
GLdB = -20 * log10(Apeak); 
freq_dbov = freq_db + GLdB;


db_max = max(max(freq_dbov));
db_range = [db_max-db_range  db_max];

imag = imagesc(t, freq, freq_dbov, db_range);
axis('xy');
xlabel('Time(s)');
ylabel('Frequency(Hz)');
colormap(flipud(gray));
%make it greyscale for clarity



function slice = getSlice (x, wind, is)
isc = max(is, 1);
siglen = length(x);
ie = min(is+length(wind)-1, siglen);
Nv = ie - isc + 1;

slice = zeros(length(wind),1);
if (Nv > 0)
    offs = isc - is;
    
    slice(offs+1:offs+Nv) = x(isc:ie);
end
slice = wind.*slice;

