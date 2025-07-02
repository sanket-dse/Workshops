%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The author of the code is Sanket S Houde working under Dr.Pragathi
% Balasubramani at Department of Cognitive Sciences, IIT Kanpur. You can
% reach me or Dr. Pragathi via our email IDs. 
% Sanket S Houde - sankethoude@gmail.com
% Dr. Pragathi Balasubramani - pbalasub@iitk.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requirements : 
% 1. Pre-condition file(.mat format)    2. Post-condition file (.mat format)
% What it does ?
% 1. Filters the EGG data 
% 2. Calculates the Power Spectral Density (PSD) using padding of the
% signals. Calculates dominant frequency and power using the PSD.
% 3. Calculates the percentage of 1 min windows which had dominant
% frequency in normogastric, tachygstric and bradygastric ranges. 
% 4. Calculates the wave propogation using the 2d fft method. 
% Reference : https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000487
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading files
%addpath(genpath('C:\Users\sanke\Downloads\NIMHANS_EGG_Workshop\adinstruments_sdk_matlab-master'));
%adi.convert('C:\Users\sanke\Downloads\NIMHANS_EGG_Workshop\HC4_post_lunch_EGG.adicht');
raw_egg_pre = load('HC3_pre lunch_EGG.mat');
raw_egg_post = load('HC4_post_lunch_EGG.mat');

%% Creating a egg signal file
egg_sig_pre = [];
egg_sig_pre(1,:) = raw_egg_pre.data__chan_1_rec_1;
egg_sig_pre(2,:) = raw_egg_pre.data__chan_2_rec_1;
egg_sig_pre(3,:) = raw_egg_pre.data__chan_3_rec_1;
egg_sig_pre(4,:) = raw_egg_pre.data__chan_4_rec_1;
egg_sig_pre(5,:) = raw_egg_pre.data__chan_5_rec_1;
egg_sig_pre(6,:) = raw_egg_pre.data__chan_6_rec_1;
egg_sig_pre(7,:) = raw_egg_pre.data__chan_7_rec_1;
egg_sig_pre(8,:) = raw_egg_pre.data__chan_8_rec_1;
egg_sig_pre = fillmissing(egg_sig_pre, 'previous'); %To take care of missing points

egg_sig_post = [];
egg_sig_post(1,:) = raw_egg_post.data__chan_1_rec_2;
egg_sig_post(2,:) = raw_egg_post.data__chan_2_rec_2;
egg_sig_post(3,:) = raw_egg_post.data__chan_3_rec_2;
egg_sig_post(4,:) = raw_egg_post.data__chan_4_rec_2;
egg_sig_post(5,:) = raw_egg_post.data__chan_5_rec_2;
egg_sig_post(6,:) = raw_egg_post.data__chan_6_rec_2;
egg_sig_post(7,:) = raw_egg_post.data__chan_7_rec_2;
egg_sig_post(8,:) = raw_egg_post.data__chan_8_rec_2;
egg_sig_post = fillmissing(egg_sig_post, 'previous'); %To take care of missing points

%% Broadband filtering the EGG channels 

srate               = 10;
center_frequency    = 0.07;        
bandwidth           = 0.08;
transition_width    = 0.15;
nyquist             = srate/2;
ffreq(1)            = 0;
ffreq(2)            = (1-transition_width)*(bandwidth - center_frequency);
ffreq(3)            = (bandwidth - center_frequency);
ffreq(4)            = (center_frequency+bandwidth);
ffreq(5)            = (1+transition_width)*(center_frequency+bandwidth);
ffreq(6)            = nyquist;
ffreq               = ffreq/nyquist;
fOrder              = 3; % in cycles
filterOrder         = fOrder*fix(srate/(bandwidth - center_frequency)); %in samples
idealresponse       = [ 0 0 1 1 0 0 ];
filterweights       = fir2(filterOrder,ffreq,idealresponse);
figure;freqz(filterweights,1,[],10);

% filter
%EGG_filt = EGG_raw;
egg_length_post = size(egg_sig_post,2);
egg_sig_padded_post = [zeros(10000,8) ; egg_sig_post' ; zeros(225100 - egg_length_post,8)];
egg_length_pre = size(egg_sig_pre,2);
egg_sig_padded_pre = [zeros(10000,8) ; egg_sig_pre' ; zeros(225100 - egg_length_pre,8)];
disp('Filtering EGG - this will take some time');
s1_filt_post   = filtfilt(filterweights,1,egg_sig_padded_post);
s1_filt_pre   = filtfilt(filterweights,1,egg_sig_padded_pre);
fprintf('filtered\n');

%Chopping off the padding and 2000 points at the start and end 
s1_filt_post = s1_filt_post(12000:egg_length_post+8000,:); 
s1_filt_pre = s1_filt_pre(12000:egg_length_pre+8000,:);

%% Plotting raw EGG and filtered EGG data

figure;
subplot(2,1,2)
plot(linspace(0,size(s1_filt_pre,1)/10.0,size(s1_filt_pre,1)),s1_filt_pre(:,4));hold on;
plot(linspace(0,size(s1_filt_post,1)/10.0,size(s1_filt_post,1)),s1_filt_post(:,4)); hold off;
ylabel('Amplitude'); xlabel('Time(s)');
xlim([0,1800])
legend({'Pre lunch','Post lunch'});
title('Filtered EGG data');

subplot(2,1,1)
plot(linspace(0,size(egg_sig_pre,2)/10.0,size(egg_sig_pre,2)),egg_sig_pre(4,:)');hold on;
plot(linspace(0,size(egg_sig_post,2)/10.0,size(egg_sig_post,2)),egg_sig_post(4,:)'); hold off;
ylabel('Amplitude'); xlabel('Time(s)');
xlim([0,2100])
legend({'Pre lunch','Post lunch'});
title('Raw EGG data')

%% Padded PSD calculation 

% This is the PSD by padding the windows for preprandial recording  
fs_egg = 10; %Hz
windows = floor(size(s1_filt_pre,1) / (60*fs_egg));
electrode = 4;
pxx_pre = [];
for i = 1:windows
    egg_sig = s1_filt_pre(((60*fs_egg)*(i-1))+1:60*fs_egg*i,electrode); %Taking a window of the signal
    egg_sig_padded = [zeros(1000,1) ; egg_sig ; zeros(1000,1)]; %Zero padding the signal
    %Calculating the psd
    fft_length = size(egg_sig_padded,1);
    [pxx_pre(:,i),freq] = pwelch(egg_sig_padded,[],[],fft_length,fs_egg);
    %figure;plot(freq(1:250),pxx_pre(1:250,i));
end

% This is the PSD by padding the windows for postprandial recording  
fs_egg = 10;
windows = floor(size(s1_filt_post,1) / (60*fs_egg));
electrode = 4; 
pxx_post = [];
for i = 1:windows
    egg_sig = s1_filt_post(((60*fs_egg)*(i-1))+1:60*fs_egg*i,electrode); %Taking a window of the signal
    egg_sig_padded = [zeros(1000,1) ; egg_sig ; zeros(1000,1)]; %Zero padding the signal
    %Calculating the psd
    fft_length = size(egg_sig_padded,1);
    [pxx_post(:,i),freq] = pwelch(egg_sig_padded,[],[],fft_length,fs_egg);
    %figure;plot(freq(1:250),pxx_post(1:250,i));
end

%% Calculating the dominant frequency and power

%Calculating the power spectrum by averaging the windows
psd_post = mean(pxx_post,2);
psd_pre = mean(pxx_pre,2);

%%Plotting the figures

% Getting the dominant power and freq
dom_power_pre = max(psd_pre);
dom_freq_pre = freq(find(psd_pre == dom_power_pre));

dom_power_post = max(psd_post);
dom_freq_post = freq(find(psd_post == dom_power_post));


figure;
plot(freq(1:50),psd_pre(1:50));hold on;
plot(freq(1:50),psd_post(1:50))
plot(dom_freq_pre,dom_power_pre,'r*');
plot(dom_freq_post,dom_power_post,'r*');hold off;
title('Power spectral density');
ylabel('Power'); xlabel('Frequency(Hz)');
legend({'Pre lunch','Post lunch'})

%{
figure;plot(freq(1:50),psd_pre(1:50));
title('Welchs periodogram of pre lunch recording ');
xlabel('Frequency (Hz)');
ylabel('Power')
txt = append('\leftarrow Dominant power = ' , num2str(dom_power));
text(dom_freq,dom_power,txt,'FontSize',14)
fprintf('Dominant Frequency in pre lunch recording = %f \n',dom_freq);

dom_power = max(psd_post);
dom_freq = freq(find(psd_post == dom_power));

figure;plot(freq(1:50),psd_post(1:50));
title('Welchs periodogram of post lunch recording ');
xlabel('Frequency (Hz)');
ylabel('Power')
txt = append('\leftarrow Dominant power = ' , num2str(dom_power));
text(dom_freq,dom_power,txt,'FontSize',14)
fprintf('Dominant Frequency in post lunch recording = %f \n ',dom_freq);
%}


%% Calculating percentage of normogastria, bradygastria, tachygastria

idx = find(freq > 0.0083 & freq < 1.5); %Finding indices in the range of interest

% For pre lunch
normo_count = 0;
tachy_count = 0;
brady_count = 0;
dom_freq = [];
dom_power = []; 

for i = 1:size(pxx_pre,2)

    dom_power(i) = max(pxx_pre(idx,i));
    dom_freq(i) = freq(find(pxx_pre(:,i) == dom_power(i)));

    if dom_freq(i) > 0.0083 & dom_freq(i) < 0.03
        brady_count = brady_count + 1;
    elseif dom_freq(i) > 0.03 & dom_freq(i) < 0.07
        normo_count = normo_count + 1;
    elseif dom_freq(i) > 0.07 & dom_freq(i) < 0.15
        tachy_count = tachy_count + 1;
    end

end

tachy_perc_pre = (tachy_count / size(pxx_pre,2))*100;
brady_perc_pre = (brady_count / size(pxx_pre,2))*100;
normo_perc_pre = (normo_count / size(pxx_pre,2))*100;

%fprintf('For pre lunch recording \n\n');
%fprintf('Tachygastric percentage = %f \n',tachy_perc);
%fprintf('Bradygastric percentage = %f \n',brady_perc);
%fprintf('Normogastric percentage = %f \n\n',normo_perc);

% For post lunch 

normo_count = 0;
tachy_count = 0;
brady_count = 0;
dom_freq = [];
dom_power = [];

for i = 1:size(pxx_post,2)

    dom_power(i) = max(pxx_post(idx,i));
    dom_freq(i) = freq(find(pxx_post(:,i) == dom_power(i)));

    if dom_freq(i) > 0.0083 & dom_freq(i) < 0.03
        brady_count = brady_count + 1;
    elseif dom_freq(i) > 0.03 & dom_freq(i) < 0.07
        normo_count = normo_count + 1;
    elseif dom_freq(i) > 0.07 & dom_freq(i) < 0.15
        tachy_count = tachy_count + 1;
    end

end

tachy_perc_post = (tachy_count / size(pxx_pre,2))*100;
brady_perc_post = (brady_count / size(pxx_pre,2))*100;
normo_perc_post = (normo_count / size(pxx_pre,2))*100;

%fprintf('For post lunch recording \n\n');
%fprintf('Tachygastric percentage = %f \n',tachy_perc);
%fprintf('Bradygastric percentage = %f \n',brady_perc);
%fprintf('Normogastric percentage = %f \n',normo_perc);

%%% Barplots
x = {'Bradygastric','Normogastric','Tachygastric'};
vals = [brady_perc_pre, brady_perc_post; normo_perc_pre, normo_perc_post; tachy_perc_pre, tachy_perc_post];
figure; bar(x,vals);
legend({'Pre lunch','Post lunch'});
title('Percentage of Normogastria/Tachygastria/Bradygastria ');
ylabel('Percentage');

%% Slow wave propagation

electrodes = [3,5,4,1]; %electrodes' order
fs_egg = 10; %Hz

%Filtering in the normogastric range 
egg_length_post = size(egg_sig_post,2);
egg_sig_padded_post = [zeros(10000,8) ; egg_sig_post' ; zeros(225100 - egg_length_post,8)];
egg_length_pre = size(egg_sig_pre,2);
egg_sig_padded_pre = [zeros(10000,8) ; egg_sig_pre' ; zeros(225100 - egg_length_pre,8)];
disp('Filtering EGG - this will take some time');
s1_filt_post = bandpass(egg_sig_padded_post,[0.03 0.07],fs_egg);
s1_filt_pre = bandpass(egg_sig_padded_pre,[0.03 0.07],fs_egg);
fprintf('filtered\n');

%Chopping off the padding and 2000 points at the start and end 
s1_filt_post = s1_filt_post(12000:egg_length_post+8000,:); 
s1_filt_pre = s1_filt_pre(12000:egg_length_pre+8000,:);

% For pre and post prandial
[FW_pre, BW_pre] = logratio_EGG(s1_filt_pre(:,electrodes)','Pre-prandial EGG data');
[FW_post, BW_post] = logratio_EGG(s1_filt_post(:,electrodes)','Post-prandial EGG data');

% Calculating Log ratios
logratio_pre = log(FW_pre./BW_pre);
logratio_post = log(FW_post./BW_post);

%Printing the log-ratios
fprintf('Log ratio for Pre-prandial : %f\n', logratio_pre);
fprintf('Log ratio for post-prandial : %f\n', logratio_post);

%% Function to find boundaries

function [thX1,thX2]=findingXboundaries(twod_fft,fx)

    %this is just a makiavelic way to find the boundaries of the interesting frequencies (i.e. thX1,thX2) ##to improve
    sumFFT=sum(twod_fft);
    threshold=0.05*max(sumFFT);
    [aa,var]=sort(abs(sumFFT-threshold));
    var=var(1:40); %40 big enough not to cut. ##to improve
    thX1=fx(min(var));
    thX2=fx(max(var));

end

%% Function to find logratios in EGG data 

function [FW, BW] = logratio_EGG(egg_sig, title_str)
   % Calculating number of 1 min epochs 
   num_epochs = size(egg_sig,2) / (60*10); % Divided by 60s*10Hz (10 Hz is the sampling rate)
   
   % Setup plot
   time_points = linspace(0,60,600); %seconds
   num_electrodes = size(egg_sig,1);
   clims = [0 5];
   
   figure(1), clf
   subplot(211)
   imageh = imagesc(rand(num_electrodes,length(time_points)));
   ylabel('Layers'); xlabel('Time (s)');
   title('Space Domain')

   subplot(223)
   amph = imagesc(rand(num_electrodes,length(time_points)));
   axis xy
   clim(clims);
   title('Amplitude spectrum')

   subplot(224)
   x = {'Forward', 'Backward'};
   y = [1 1];
   phaseh = bar(x,y);
   axis xy
   title('Propagation strengths')
   
   tempFW = [];
   tempBW = [];

   for epoch = 1:num_epochs
                 
       % Computing 2d fft
       egg_sig_1min = egg_sig(:,1+(epoch-1)*60*10:epoch*60*10);
       twod_fft = abs(fftshift(fft2(egg_sig_1min)));

       a=1; %width pixel to compute the frequency

       [m,n]=size(egg_sig_1min); %numberPixel x and y
       
       samplingRate = 10; %Hertz
       durationSignal=n/samplingRate; %it's in second.
       dF = 1/durationSignal;
       fx = -samplingRate/2:dF:samplingRate/2-dF;
       fy=((1)/2)*linspace(-1,1,m);

       [thX1,thX2]=findingXboundaries(twod_fft,fx);

       upperPart=twod_fft(fy>0,(fx<0 & fx>thX1));
       lowerPart=twod_fft(fy<0,(fx<0 & fx>thX1));

       tempFW(epoch) = mean(upperPart(:));
       tempBW(epoch) = mean(lowerPart(:));

       %update plots
       set(imageh, 'CData', egg_sig_1min);
       set(amph, 'CData', twod_fft);
       set(phaseh, 'YData', [tempFW(epoch), tempBW(epoch)]);
       sgtitle(title_str)
       pause(1.0)

   end

       FW = mean(tempFW);
       BW = mean(tempBW);
end
    
%% Reserve code 
%{
%% dPLI between two electrodes

% Filter signal in the normogastric range

fprintf('Being filtered\n')
egg_length_post = size(egg_sig_post,2);
egg_sig_padded_post = [zeros(10000,8) ; egg_sig_post' ; zeros(10000,8)];
s1_filt_post = bandpass(egg_sig_padded_post,[0.03 0.07],10);
%Chopping off the padding and 2000 points at the end and the start
s1_filt_post = s1_filt_post(12000:egg_length_post+8000,:);

egg_length_pre = size(egg_sig_pre,2);
egg_sig_padded_pre = [zeros(10000,8) ; egg_sig_pre' ; zeros(10000,8)];
s1_filt_pre = bandpass(egg_sig_padded_pre,[0.03 0.07],10);
%Chopping off the padding and 2000 points at the end and the start
s1_filt_pre = s1_filt_pre(12000:egg_length_pre+8000,:);

fprintf('Filtering done\n')

%%% Calculating dPLI in the pre lunch recording
e1 = 2; %Electrode
e2 = 6;
fs_egg = 10;
windows = floor(size(s1_filt_pre,1) / (60*fs_egg));
dPLI_pre = [];

for j = 1:windows
    egg_sig = s1_filt_pre(((60*fs_egg)*(j-1))+1:60*fs_egg*j,[e1,e2]); %Taking a window of the signal
    
    % Hilbert transform
    eggH = hilbert(egg_sig);
    angle_egg = angle(eggH);

    % directed PLI
    phase_diff_TS = imag(exp(1i*(angle_egg(:,1)-angle_egg(:,2))));
    dPLI_pre(j) = mean(heaviside(phase_diff_TS));
end
%parameters
back_perc_pre = (sum(dPLI_pre < 0.5)/windows)*100;
forward_perc_pre = (sum(dPLI_pre > 0.5)/windows)*100;
mean_forward_pre = mean(dPLI_pre(dPLI_pre > 0.5));
mean_back_pre = mean(dPLI_pre(dPLI_pre < 0.5));

%%% Calculating dPLI in the post lunch recording
fs_egg = 10;
windows = floor(size(s1_filt_post,1) / (60*fs_egg));
dPLI_post = [];

for j = 1:windows
    egg_sig = s1_filt_post(((60*fs_egg)*(j-1))+1:60*fs_egg*j,[e1,e2]); %Taking a window of the signal
    
    % Hilbert transform
    eggH = hilbert(egg_sig);
    angle_egg = angle(eggH);

    % directed PLI
    phase_diff_TS = imag(exp(1i*(angle_egg(:,1)-angle_egg(:,2))));
    dPLI_post(j) = mean(heaviside(phase_diff_TS));
end

%parameters
back_perc_post = (sum(dPLI_post < 0.5)/windows)*100;
forward_perc_post = (sum(dPLI_post > 0.5)/windows)*100;
mean_forward_post = mean(dPLI_post(dPLI_post > 0.5));
mean_back_post = mean(dPLI_post(dPLI_post < 0.5));

%%% Barplots for showing the metrics

x = {'Back propagation','Forward propagation'};
vals = [back_perc_pre, back_perc_post; forward_perc_pre,forward_perc_post];
figure; bar(x,vals);
legend({'Pre lunch','Post lunch'});
title('Normogastric wave propogation percentage');
ylabel('Propagation percentage');

x = {'Back propagation','Forward propagation'};
vals = [mean_back_pre, mean_back_post; mean_forward_pre,mean_forward_post];
figure; bar(x,vals);hold on;
yline(0.5,'--');
legend({'Pre lunch','Post lunch'});
title('Strength of propogation');
ylabel('dPLI');

%% Slow wave propogation

fs_egg = 10;
e1 = 2;
e2 = 6;

%Post lunch recording
windows = floor(size(s1_filt_post,1) / (60*fs_egg));
lag_win_count_post = 0;
lag_strength_post = [];
lead_win_count_post = 0;
lead_strength_post = [];
lag_length = 3*fs_egg;


for i=1:windows
    [c,lags] = xcov(s1_filt_post(((60*fs_egg)*(i-1))+1:60*fs_egg*i,e1),s1_filt_post(((60*fs_egg)*(i-1))+1:60*fs_egg*i,e2));
    %figure;stem(lags,c);
    if(lags(find(c==max(c))) > lag_length)
        lead_win_count_post = lead_win_count_post + 1;
        lead_strength_post(lead_win_count_post) = abs(lags(c==max(c)));
    end    
    
    if(lags(find(c==max(c))) < -lag_length)
        lag_win_count_post = lag_win_count_post + 1;
        lag_strength_post(lag_win_count_post) = abs(lags(c==max(c)));
    end
    
end
lag_perc_post = lag_win_count_post / windows;
lead_perc_post = lead_win_count_post / windows;

%Pre lunch recording
fs_egg = 10;
windows = floor(size(s1_filt_pre,1) / (60*fs_egg));
lag_win_count_pre = 0;
lag_strength_pre = [];
lead_win_count_pre = 0;
lag_strength_pre = [];

for i=1:windows
    [c,lags] = xcov(s1_filt_pre(((60*fs_egg)*(i-1))+1:60*fs_egg*i,e1),s1_filt_pre(((60*fs_egg)*(i-1))+1:60*fs_egg*i,e2));
    %figure;stem(lags,c);
    if(lags(find(c==max(c))) > lag_length)
        lead_win_count_pre = lead_win_count_pre + 1;
        lead_strength_pre(lead_win_count_pre) = abs(lags(c==max(c)));
    end    
    
    if(lags(find(c==max(c))) < -lag_length)
        lag_win_count_pre = lag_win_count_pre + 1;
        lag_strength_pre(lag_win_count_pre) = abs(lags(c==max(c)));
    end
end
lag_perc_pre = lag_win_count_pre / windows;
lead_perc_pre = lead_win_count_pre / windows;

% barplots 
x = {'Back propagation','Forward propagation'};
vals = [lag_perc_pre, lag_perc_post; lead_perc_pre,lead_perc_post];
figure; bar(x,vals);
legend({'Pre lunch','Post lunch'});
title('Normogastric wave propogation percentage');
ylabel('Propagation percentage');

x = {'Back propagation','Forward propagation'};
vals = [mean(lag_strength_pre)/fs_egg, mean(lag_strength_post)/fs_egg; mean(lead_strength_pre)/fs_egg,mean(lead_strength_post)/fs_egg];
figure; bar(x,vals);
legend({'Pre lunch','Post lunch'});
title('Strength of propogation');
ylabel('Averaged time lag(s)');
%}