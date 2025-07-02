%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The author of the code is Sanket S Houde working under Dr.Pragathi
% Balasubramani at Department of Cognitive Sciences, IIT Kanpur. You can
% reach me or Dr. Pragathi via our email IDs. 
% Sanket S Houde - sankethoude@gmail.com
% Dr. Pragathi Balasubramani - pbalasub@iitk.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This creates two synthetic signals : high frequency representing EEG waves and low frequency
% waves representing EGG waves. Then the phase amplitude coupling is calculated between the two 
% waves. You are better off using the applet. This is to give an inner
% working of the app.
% Reference : https://www.sciencedirect.com/science/article/pii/S0165027011004730?via%3Dihub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adding important folders to path 
addpath(genpath('C:\Users\sanke\Downloads\PD_EEG_EGG\PACT-master'));

%% Creating EEG and EGG synthetic waves 

%%Time specifications:
Fs = 1000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 3;             % seconds
t = (0:dt:StopTime-dt)';     % seconds


%%Slow wave

% parameters : 1) intermediate PAC : vertical_shift = 0.5, horizontal_shift = pi/3, 
% Amp = 0.5, Second slow_waves_egg expression
% 2) High PAC : vertical_shift = 0 , horizontal_shift = pi/2 , Amp=1, first slow_waves_egg
% expression
% 3) Low PAC : vertical_shift = 0, horizontal_shift = pi, Amp = 0.5 , second slow_waves_egg
% expression

%%%%%% Creation of slow EGG waves
Fc = 1.5; % hertz
Amp = 1;
horizontal_shift = pi/2;
vertical_shift = 0;
slow_waves_egg = abs(Amp*cos(2*pi*Fc*t + horizontal_shift))- vertical_shift; %first slow_waves_egg expression
%slow_waves_egg = Amp*cos(2*pi*Fc*t + horizontal_shift) + vertical_shift; %second slow_waves_egg expression

slow_sim = (slow_waves_egg);
%%%%%%%

%%%%%%%%% Creation of fast EEG waves

% Slow wave for mixing
Fc = 1.5;  % hertz
Amp = 1;
phaseAng = pi/2;
slow_waves = Amp*cos(2*pi*Fc*t + phaseAng);

% Fast wave:
Fc = 30;  % hertz
Amp = 1;
phaseAng = pi/2;
high_waves = Amp*cos(2*pi*Fc*t + phaseAng);

fast_sim = (slow_waves.*high_waves);
%%%%%%%%%%%%



%% Calculating PAC

eegH = hilbert(fast_sim);
eggH = hilbert(slow_sim);

eeg_amp = abs(eegH);
egg_phase = angle(eggH);
eeg_egg_coupled_signal = eeg_amp.*exp(1i*egg_phase);

percentile_threshold = prctile(eeg_amp, 95);
indices_top_5_percentile = find(eeg_amp >= percentile_threshold);

%Overall PAC
pac_numerator = abs(sum(eeg_egg_coupled_signal));
pac_denominator = sqrt(sum(eeg_amp.*eeg_amp))*sqrt(length(high_waves));
pac_overall = pac_numerator/pac_denominator;

%Top 5%ile PAC
pac_numerator = abs(sum(eeg_egg_coupled_signal(indices_top_5_percentile)));
pac_denominator = sqrt(sum(eeg_amp(indices_top_5_percentile).*eeg_amp(indices_top_5_percentile)))*sqrt(length(high_waves(indices_top_5_percentile)));
pac_top5percentile = pac_numerator/pac_denominator;


%% Creating the plots for PAC

phase = angle(eeg_egg_coupled_signal);
magnitude = abs(eeg_egg_coupled_signal);
angle_PAC = angle(sum(eeg_egg_coupled_signal));
denom = sqrt((eeg_amp.*eeg_amp));
magnitude = abs(eeg_egg_coupled_signal)./denom;

% Plotting the time domain signals
figure;
subplot(121)
plot(t,fast_sim);hold on;
plot(t,slow_sim,'r');hold off;
legend({'Simulated fast waves','Simulated slow waves'})

% Creating a polar plot
subplot(222)
nbins = 49;
circ_plot(phase,'hist',[], nbins,false,true); %'linewidth',0.001)

% Creating figure for power envelope of fast signal and time domain signal
% of slow signal 
subplot(224)
plot(slow_sim,'b--'); hold on;
plot(eeg_amp,'r');
legend({'Slow wave','Power envelope of fast wave'})
sgtitle(append('PAC : ',num2str(pac_overall)))

%% Reserve

%figure;plot(egg_phase,'b');hold on;
%plot(eeg_amp,'r');hold off;

%figure;
%polarplot([repelem(0,length(phase)); phase'], ...
%    [repelem(0,length(phase));magnitude'],'b'); hold on;
%polarplot([0;angle_PAC],[0;pac_overall],'r','LineWidth',3);
%legend({'Vectors',repmat({''},length(phase)-1),'Mean Vector Length'})
