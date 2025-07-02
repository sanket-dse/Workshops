%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The author of the code is Sanket S Houde working under Dr.Pragathi
% Balasubramani at Department of Cognitive Sciences, IIT Kanpur. You can
% reach me or Dr. Pragathi via our email IDs. 
% Sanket S Houde - sankethoude@gmail.com
% Dr. Pragathi Balasubramani - pbalasub@iitk.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adding important folders to path 

addpath(genpath('C:\Users\sanke\Downloads\ACCS_Workshop_2023\eeglab2023.1'));
addpath('C:\Users\sanke\Downloads\ACCS_Workshop_2023\Tasks');
addpath('C:\Users\sanke\Downloads\ACCS_Workshop_2023\xdfimport1.13');
load('Task1.mat');
%load('Task2.mat')

%% Selecting the egg electrode for gastric wave

s1_filt = bandpass(egg_data',[0.03,0.07],250); %Since the gastric slow wave is present in that freq

% Calculate power
fs_egg = 250;
nfft = 60*fs_egg;
window = hanning(nfft); %rectwin()
noverlap = nfft/2;
[Pxx1,w1] =  pwelch(s1_filt,window,noverlap,nfft,fs_egg);
power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

range_idx = find(w1>0.03 & w1<0.07);
mean_power = sum(power_norm(range_idx,:))./sum(power_norm);
egg_electrode = find(mean_power == max(mean_power)); 

egg_sig = s1_filt(:,egg_electrode); %selecting the EGG electrode with the highest power

%% Visualizing the EGG electrodes power spectrum

% Calculate power
fs_egg = 250;
nfft = 60*fs_egg;
window = hanning(nfft); %rectwin()
noverlap = nfft/2;
[Pxx1,w1] =  pwelch(egg_data',window,noverlap,nfft,fs_egg);
power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

figure;
for electrode = 1:size(egg_data,1)
    plot(w1(1:120),log(power_norm(1:120,electrode)));
    hold on;
end
legend({'5','6','7','8'})
title('Power spectrum for EGG')
%% Selecting the heart signal electrode
heart_filt = bandpass(egg_data',[1.2,1.8],250); % Put the range of cardiac rhythm into this line

% Calculate power
fs_egg = 250;
nfft = 60*fs_egg;
window = hanning(nfft); %rectwin()
noverlap = nfft/2;
[Pxx1,w1] =  pwelch(heart_filt,window,noverlap,nfft,fs_egg);
power_norm = Pxx1./ sum(Pxx1); %normalizing the PSD 

range_idx = find(w1>1.2 & w1<1.6);
mean_power = sum(power_norm(range_idx,:))./sum(power_norm);
heart_electrode = find(mean_power == max(mean_power)); 

heart_sig = heart_filt(:,heart_electrode);

%% preprocessing EEG data
EEG_sanket = load_xdf('C:\Users\sanke\Downloads\ACCS_Workshop_2023\1001 SH_T1.xdf');

%%% getting channel locations (old) : first run this 

nchans = 32;
chanlocs = [];

for i=1:nchans 

    
    chanlocs(:,i).X = str2double(EEG_sanket{1, 3}.info.desc.channels.channel{1, i}.location.X);
    chanlocs(:,i).Y = str2double(EEG_sanket{1, 3}.info.desc.channels.channel{1, i}.location.Z);
    chanlocs(:,i).Z = str2double(EEG_sanket{1, 3}.info.desc.channels.channel{1, i}.location.Y);
    
    
    %Use cart2topo = cart2sph() -> sph2topo().
    % First convert cartesian to spherical coordinates and then to topoplot
    % coordinates

    [sph_theta, sph_phi, chanlocs(:,i).sph_radius] = cart2sph(chanlocs(:,i).X,chanlocs(:,i).Y,chanlocs(:,i).Z);
    chanlocs(:,i).sph_phi = rad2deg(sph_phi);
    chanlocs(:,i).sph_theta = rad2deg(sph_theta);

    [chanlocs(:,i).urchanlocs] = i;
    [~,chanlocs(:,i).theta,chanlocs(:,i).radius] = sph2topo([chanlocs(:,i).urchanlocs,chanlocs(:,i).sph_phi,chanlocs(:,i).sph_theta]);
    chanlocs(:,i).labels = EEG_sanket{1, 3}.info.desc.channels.channel{1, i}.label;

end

%%% getting channel locations (new)

%reading the chanlocs txt file 

chansinfo = readmatrix('C:\Users\sanke\Downloads\ACCS_Workshop_2023\chanlocs.txt');
temp = readtable('C:\Users\sanke\Downloads\ACCS_Workshop_2023\chanlocs.txt');
chanlabels = table2array(temp(:,1));
chansinfo = chansinfo(:,2:end);
chanlocs_new = [];
nchans = 32;

for i=1:nchans

    idx = find(strcmpi(chanlabels,chanlocs(:,i).labels));
    chanlocs_new(:,i).X = chansinfo(idx,1);
    chanlocs_new(:,i).Y = chansinfo(idx,2);
    chanlocs_new(:,i).Z = chansinfo(idx,3);

    [sph_theta, sph_phi, chanlocs_new(:,i).sph_radius] = cart2sph(chanlocs_new(:,i).X,chanlocs_new(:,i).Y,chanlocs_new(:,i).Z);
    chanlocs_new(:,i).sph_phi = rad2deg(sph_phi);
    chanlocs_new(:,i).sph_theta = rad2deg(sph_theta);

    [chanlocs_new(:,i).urchanlocs] = i;
    [~,chanlocs_new(:,i).theta,chanlocs_new(:,i).radius] = sph2topo([chanlocs_new(:,i).urchanlocs,chanlocs_new(:,i).sph_phi,chanlocs_new(:,i).sph_theta]);
    chanlocs_new(:,i).labels = char(chanlabels(idx));

    chanlocs_new(:,i).theta = chanlocs_new(:,i).theta + 90 ; %this was done because the positions were rotated by 90 degrees


end

%%% Creating the EEG struct

EEG_raw = pop_importdata('dataformat','matlab','nbchan',0,'data',eeg_data, ...
        'srate',250,'pnts',0,'xmin',0);

EEG_raw.trials = 1;
EEG_raw.nbchan = size(EEG_raw.data,1);
EEG_raw.pnts = size(EEG_raw.data,2);
EEG_raw.srate = 250;
EEG_raw.xmin = 0;
EEG_raw.xmax = size(EEG_raw.data,2)/EEG_raw.srate;
EEG_raw.times = linspace(EEG_raw.xmin,EEG_raw.xmax,EEG_raw.pnts);
EEG_raw.etc = [];

%reading channel locations 
EEG_raw.chanlocs = chanlocs_new;

%%% preprocessing data - FIR filter

EEG_filt = pop_eegfiltnew(EEG_raw, 'locutoff',  1, 'hicutoff',  45, 'filtorder', 9000, 'plotfreqz', 0); % this gives better results

%%% Remove bad channels using clean_channels function
EEG_filt_rmchan = clean_channels(EEG_filt);

%%% Performing ICA
ch_list = 1:EEG_filt_rmchan.nbchan;
EEG_filt_rmchan = pop_runica(EEG_filt_rmchan,'runica');
EEG_filt_rmchan.icachansind = double(1:EEG_filt_rmchan.nbchan);
EEG_filt_rmchan = iclabel(EEG_filt_rmchan);
    
fprintf('ICA done\n')

%%% Remove bad components and reconstruct data
threshold_signal = 0.05; % brain component with less than 5% confidence is removed
cls = EEG_filt_rmchan.etc.ic_classification.ICLabel.classes;
cls_score = EEG_filt_rmchan.etc.ic_classification.ICLabel.classifications;
bad_comp = [];
for cmp=1:size(EEG_filt_rmchan.icachansind,2)
    if cls_score(cmp,1)<threshold_signal
        bad_comp = [bad_comp,cmp];
    end
end

EEG_ica = pop_subcomp(EEG_filt_rmchan, bad_comp, 0);
    
fprintf('Bad components removed\n');

%%% Identifying bad portions of the data and reconstructing using ASR

%EEG_clean = clean_artifacts(EEG_ica);
%fprintf('Bad portions of data identified and reconstructed \n');

%%% Interpolate bad channels
EEG_interpol = pop_interp(EEG_ica, EEG_raw.chanlocs, 'spherical');
fprintf('Bad channels interpolated\n');

%%% Re-reference the data
EEG_reref = pop_reref(EEG_interpol,[]);
fprintf('Data rereferenced\n');

%% Sectioning off EGG using the event markers
%{
num_events = length(EEG_reref.event);

for event = 1:num_events
    time_pnt = EEG_reref.event(1,event).latency;
    num_pnts = EEG_reref.event(1,event).duration
    egg_sig(time_pnt:time_pnt+num_pnts-1,:) = nan;
end
egg_sig = rmmissing(egg_sig);
%}
%%  Get the median EEG data for the brain regions

%Defining electrodes for different brain regions
frontal = {'fp1','fp2','f3','f4','fz'};
central = {'fc1','fc2','c3','c4','cz'};
parietal = {'p3','p4','cp1','cp2','pz'};
occipital = {'poz','o1','o2'};
left_temporal = {'f7','fc5','t7','cp5','p7','tp9','ft9'};
right_temporal = {'f8','fc6','cp6','p8','ft10','t8','tp10'};

frontal_data=[];
occipital_data = [];
parietal_data = [];
central_data = [];
left_temporal_data = [];
right_temporal_data = [];
EEG_data = [];
            
for i=1:size(frontal,2)
    frontal_data(i,:) = EEG_reref.data(find(strcmpi({EEG_reref.chanlocs.labels},frontal{i})),:);
end

for i=1:size(occipital,2)
    occipital_data(i,:) = EEG_reref.data(find(strcmpi({EEG_reref.chanlocs.labels},occipital{i})),:);
end

for i=1:size(parietal,2)
    parietal_data(i,:) = EEG_reref.data(find(strcmpi({EEG_reref.chanlocs.labels},parietal{i})),:);
end

for i=1:size(central,2)
    central_data(i,:) = EEG_reref.data(find(strcmpi({EEG_reref.chanlocs.labels},central{i})),:);
end

for i=1:size(left_temporal,2)
    left_temporal_data(i,:) = EEG_reref.data(find(strcmpi({EEG_reref.chanlocs.labels},left_temporal{i})),:);
end

for i=1:size(right_temporal,2)
    right_temporal_data(i,:) = EEG_reref.data(find(strcmpi({EEG_reref.chanlocs.labels},right_temporal{i})),:);
end

frontal_data = median(frontal_data,1);
occipital_data = median(occipital_data,1);
parietal_data = median(parietal_data,1);
central_data = median(central_data,1);
left_temporal_data = median(left_temporal_data,1);
right_temporal_data = median(right_temporal_data,1);

EEG_data = [frontal_data;occipital_data; parietal_data; central_data; left_temporal_data; right_temporal_data];

%% Calculating and plotting power in different bands

fs = EEG_reref.srate;
nfft = 5*fs;
window = hanning(nfft); %rectwin()
noverlap = nfft/2;
[Pxx,w] =  pwelch(EEG_reref.data',window,noverlap,nfft,fs);
power_norm = Pxx./sum(Pxx);

b_ind = find(w>13 & w<30);
beta_mean_pow = mean(power_norm(b_ind,:));
beta_tot_pow = sum(power_norm(b_ind,:));
beta_relative_pow = sum(power_norm(b_ind,:))./sum(power_norm);

b_ind = find(w>4 & w<8);
theta_mean_pow = mean(power_norm(b_ind,:));

b_ind = find(w>8 & w<13);
alpha_mean_pow = mean(power_norm(b_ind,:));

figure;
topoplot(beta_mean_pow,EEG_reref.chanlocs);
colorbar;
caxis([min(beta_mean_pow),max(beta_mean_pow)])
title('Beta power - Breathe in Breathe out')

figure;
topoplot(alpha_mean_pow,EEG_reref.chanlocs);
colorbar;
caxis([min(alpha_mean_pow),max(alpha_mean_pow)])
title('Alpha power - Breathe in Breathe out')

figure;
topoplot(theta_mean_pow,EEG_reref.chanlocs);
colorbar;
caxis([min(theta_mean_pow),max(theta_mean_pow)])
title('Theta power - Breathe in Breathe out task')

%%  Calculating PAC between alpha, beta and theta bands of EEG with EGG

EEG_beta = bandpass(EEG_data',[13,30],250);
EEG_alpha = bandpass(EEG_data',[8,13],250);
EEG_theta = bandpass(EEG_data',[4,8],250);

EGG_norm = (egg_sig-mean(egg_sig,"all"))/std(egg_sig,[],"all"); %Normalizing the EGG data

pac_gastric = [];
pac_gastric_top5 = [];

for freq_band = 1:3
    if freq_band==1
        EEG = EEG_theta;
    elseif freq_band==2
        EEG = EEG_alpha;
    else
        EEG = EEG_beta;
    end

    for region = 1:6
        % get coupled signal
        pac_amp = abs(hilbert(EEG(:,region)));
        egg_phase = angle(hilbert(EGG_norm));
        egg_eeg_coupled_signal = pac_amp.*exp(1i*egg_phase);
    
        % get top 5 percentile indices of amplitude
        percentile_threshold = prctile(pac_amp, 95);
        indices_top_5_percentile = find(pac_amp >= percentile_threshold);
    
        % compute overall PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal));
        pac_denominator = sqrt(sum(pac_amp.*pac_amp))*sqrt(size(EEG_beta,1));
        pac_gastric(freq_band,region) = pac_numerator/pac_denominator;
    
        % compute top 5 percentile PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal(indices_top_5_percentile)));
        pac_denominator = sqrt(sum(pac_amp(indices_top_5_percentile).*pac_amp(indices_top_5_percentile)))*sqrt(size(EEG_beta(indices_top_5_percentile),1));
        pac_gastric_top5(freq_band,region) = pac_numerator/pac_denominator;

        indices_top5(:,freq_band,region) = indices_top_5_percentile; % Storing the top5 percentile amplitude indices 
    end
end

%% Calculating PAC between alpha, beta and theta bands of EEG with heart signal

EEG_beta = bandpass(EEG_data',[13,30],250);
EEG_alpha = bandpass(EEG_data',[8,13],250);
EEG_theta = bandpass(EEG_data',[4,8],250);

heart_norm = (heart_sig-mean(heart_sig,"all"))/std(heart_sig,[],"all"); %Normalizing the EGG data

pac_heart = [];
pac_heart_top5 = [];

for freq_band = 1:3
    if freq_band==1
        EEG = EEG_theta;
    elseif freq_band==2
        EEG = EEG_alpha;
    else
        EEG = EEG_beta;
    end

    for region = 1:6
        % get coupled signal
        pac_amp = abs(hilbert(EEG(:,region)));
        egg_phase = angle(hilbert(heart_norm));
        egg_eeg_coupled_signal = pac_amp.*exp(1i*egg_phase);
    
        % get top 5 percentile indices of amplitude
        percentile_threshold = prctile(pac_amp, 95);
        indices_top_5_percentile = find(pac_amp >= percentile_threshold);
    
        % compute overall PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal));
        pac_denominator = sqrt(sum(pac_amp.*pac_amp))*sqrt(size(EEG_beta,1));
        pac_heart(freq_band,region) = pac_numerator/pac_denominator;
    
        % compute top 5 percentile PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal(indices_top_5_percentile)));
        pac_denominator = sqrt(sum(pac_amp(indices_top_5_percentile).*pac_amp(indices_top_5_percentile)))*sqrt(size(EEG_beta(indices_top_5_percentile),1));
        pac_heart_top5(freq_band,region) = pac_numerator/pac_denominator;
    end
end

%% Saving data 
%save('PAC_values/PAC_T2.mat','pac_heart_top5','pac_gastric_top5');

%% Visualizing EEG and EGG coupled analytical signal

% Just select the one with the highest top5_gastric_PAC 

figure;
idx = indices_top5(:,2,6);
pac_amp = abs(hilbert(EEG_alpha(:,6)));
egg_phase = angle(hilbert(EGG_norm));
eeg_egg_coupled_signal = pac_amp.*exp(1i*egg_phase);
phase = angle(eeg_egg_coupled_signal(idx));
denom = sqrt((pac_amp.*pac_amp));
magnitude = abs(eeg_egg_coupled_signal(idx))./denom(idx); % magnitude of all the coupled signal vectors
angle_PAC = angle(sum(eeg_egg_coupled_signal));


polarplot([repelem(0,length(phase)); phase'], ...
    [repelem(0,length(phase));magnitude'],'b'); hold on;
polarplot([0;angle_PAC],[0;pac_gastric_top5(2,6)],'r','LineWidth',3);

%% Visualizing power envelope of EEG and EGG signal

figure;
plot(indices_top5(:,2,6),pac_amp(indices_top5(:,2,6)));hold on;
plot(EGG_norm);hold on;
plot(pac_amp);

%%  Calculating PAC between alpha, beta and theta bands of EEG with EGG (For topoplot)

EEG_beta = bandpass(EEG_reref.data',[13,30],250);
EEG_alpha = bandpass(EEG_reref.data',[8,13],250);
EEG_theta = bandpass(EEG_reref.data',[4,8],250);

EGG_norm = (egg_sig-mean(egg_sig,"all"))/std(egg_sig,[],"all"); %Normalizing the EGG data

pac = [];
pac_top5 = [];

for freq_band = 1:3
    if freq_band==1
        EEG = EEG_theta;
    elseif freq_band==2
        EEG = EEG_alpha;
    else
        EEG = EEG_beta;
    end

    for region = 1:32
        % get coupled signal
        pac_amp = abs(hilbert(EEG(:,region)));
        egg_phase = angle(hilbert(EGG_norm));
        egg_eeg_coupled_signal = pac_amp.*exp(1i*egg_phase);
    
        % get top 5 percentile indices of amplitude
        percentile_threshold = prctile(pac_amp, 95);
        indices_top_5_percentile = find(pac_amp >= percentile_threshold);
    
        % compute overall PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal));
        pac_denominator = sqrt(sum(pac_amp.*pac_amp))*sqrt(size(EEG_beta,1));
        pac(freq_band,region) = pac_numerator/pac_denominator;
    
        % compute top 5 percentile PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal(indices_top_5_percentile)));
        pac_denominator = sqrt(sum(pac_amp(indices_top_5_percentile).*pac_amp(indices_top_5_percentile)))*sqrt(size(EEG_beta(indices_top_5_percentile),1));
        pac_top5(freq_band,region) = pac_numerator/pac_denominator;

        %indices_top5(:,freq_band,region) = indices_top_5_percentile; % Storing the top5 percentile amplitude indices 
    end
end

figure;
topoplot(pac_top5(3,:),EEG_reref.chanlocs);
colorbar;
caxis([min(pac_top5(3,:)),max(pac_top5(3,:))])
title('Beta gastric PAC - Breathe in breathe out')

figure;
topoplot(pac_top5(2,:),EEG_reref.chanlocs);
colorbar;
caxis([min(pac_top5(2,:)),max(pac_top5(2,:))])
title('Alpha gastric PAC - Breathe in breathe out')

figure;
topoplot(pac_top5(1,:),EEG_reref.chanlocs);
colorbar;
caxis([min(pac_top5(1,:)),max(pac_top5(1,:))])
title('Theta gastric PAC - Breathe in breathe out')

%% Calculating PAC between alpha, beta and theta bands of EEG with heart signal (For topoplot)

EEG_beta = bandpass(EEG_reref.data',[13,30],250);
EEG_alpha = bandpass(EEG_reref.data',[8,13],250);
EEG_theta = bandpass(EEG_reref.data',[4,8],250);

heart_norm = (heart_sig-mean(heart_sig,"all"))/std(heart_sig,[],"all"); %Normalizing the EGG data

pac_heart = [];
pac_heart_top5 = [];

for freq_band = 1:3
    if freq_band==1
        EEG = EEG_theta;
    elseif freq_band==2
        EEG = EEG_alpha;
    else
        EEG = EEG_beta;
    end

    for region = 1:32
        % get coupled signal
        pac_amp = abs(hilbert(EEG(:,region)));
        egg_phase = angle(hilbert(heart_norm));
        egg_eeg_coupled_signal = pac_amp.*exp(1i*egg_phase);
    
        % get top 5 percentile indices of amplitude
        percentile_threshold = prctile(pac_amp, 95);
        indices_top_5_percentile = find(pac_amp >= percentile_threshold);
    
        % compute overall PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal));
        pac_denominator = sqrt(sum(pac_amp.*pac_amp))*sqrt(size(EEG_beta,1));
        pac(freq_band,region) = pac_numerator/pac_denominator;
    
        % compute top 5 percentile PAC
        pac_numerator = abs(sum(egg_eeg_coupled_signal(indices_top_5_percentile)));
        pac_denominator = sqrt(sum(pac_amp(indices_top_5_percentile).*pac_amp(indices_top_5_percentile)))*sqrt(size(EEG_beta(indices_top_5_percentile),1));
        pac_top5(freq_band,region) = pac_numerator/pac_denominator;
    end
end


figure;
topoplot(pac_top5(3,:),EEG_reref.chanlocs);
colorbar;
caxis([min(pac_top5(3,:)),max(pac_top5(3,:))])
title('Beta heart PAC - Breathe in breathe out')

figure;
topoplot(pac_top5(2,:),EEG_reref.chanlocs);
colorbar;
caxis([min(pac_top5(2,:)),max(pac_top5(2,:))])
title('Alpha heart PAC - Breathe in breathe out')

figure;
topoplot(pac_top5(1,:),EEG_reref.chanlocs);
colorbar;
caxis([min(pac_top5(1,:)),max(pac_top5(1,:))])
title('Theta heart PAC - Breathe in breathe out')

