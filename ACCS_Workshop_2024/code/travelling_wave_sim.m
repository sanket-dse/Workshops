%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The author of the code is Sanket S Houde working under Dr.Pragathi
% Balasubramani at Department of Cognitive Sciences, IIT Kanpur. You can
% reach me or Dr. Pragathi via our email IDs. 
% Sanket S Houde - sankethoude@gmail.com
% Dr. Pragathi Balasubramani - pbalasub@iitk.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This showcases simulation of travelling waves using synthetic signals.
% There are 2 kinds of disruption : 
% 1. Layer 2. Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

sinephase = linspace(-pi/2,pi/2,50);
sinefreq = 0.05; %Hz (arbitrary units)
amplitude = 100;

%sine wave initialization
time_points = linspace(0,100,100); %seconds
num_electrodes = 100;
clims = [0 5000];

% Setup plot
figure(1), clf
subplot(121)
imageh = imagesc(rand(length(time_points)));
ylabel('Layers'); xlabel('Time (s)');
axis square, axis xy
title('Space Domain')

subplot(222)
amph = imagesc(rand(length(time_points)));
axis xy
%set(gca, 'xlim', [lims(2)-30 lims(2)+30], 'ylim', [lims(2)-30 lims(2)+30],'clim',clim);
clim(clims);
title('Amplitude spectrum')

subplot(224)
x = {'Forward', 'Backward'};
y = [1 1];
phaseh = bar(x,y);
axis xy
%axis off, axis xy
%set(gca, 'xlim', [lims(2)-30, lims(2)+30], 'ylim', [lims(2)-30, lims(2)+30])
title('Propagation strengths')

%% Loop to go through all the phases in case of a normal travelling wave

for si = 1:length(sinephase)
    
    %Creating the 2d sine wave
    %xp = x*cos(sinephase(si)) + y*sin(sinephase(si));
    %img = sin(2*pi*sinefreq*xp);
    
    %Creating the 2d sine wave
    img = [];
    for ei = 1:num_electrodes
        img(ei,:) = amplitude*sin(2*pi*sinefreq*time_points + (ei-1)*sinephase(si));
    end
    
    %Getting the 2d fft
    imgx = fftshift(fft2(img));
    powr2 = abs(imgx);
    %phas2 = angle(imgx);
    [FW, BW] = propagation_fun(img);

    %update plots
    set(imageh, 'CData', img);
    set(amph, 'CData', powr2);
    set(phaseh, 'YData', [FW, BW]);
    sgtitle('Normal Propagation')
    pause(1.0)
end

%
%% Loop to go through in case of #1 kind of disruption

for si = 1:length(sinephase)
    
    %Creating the 2d sine wave
    %xp = x*cos(sinephase(si)) + y*sin(sinephase(si));
    %img = sin(2*pi*sinefreq*xp);
    
    %Creating the 2d sine wave
    img = [];
    for ei = 1:num_electrodes/4
        img(ei,:) = amplitude*sin(2*pi*sinefreq*time_points + (ei-1)*sinephase(si));
    end

    for ei = (num_electrodes/4)+1:num_electrodes
        img(ei,:) = amplitude*sin(2*pi*sinefreq*time_points + -1*(ei-1)*sinephase(si));
    end

    %Getting the 2d fft
    imgx = fftshift(fft2(img));
    powr2 = abs(imgx);
    %phas2 = angle(imgx);
    [FW, BW] = propagation_fun(img);

    %update plots
    set(imageh, 'CData', img);
    set(amph, 'CData', powr2);
    set(phaseh, 'YData', [FW, BW]);
    sgtitle('#1 kind of disruption')
    pause(1.0)
end

%% Loop to go through in case of #2 kind of disruption

for si = 1:length(sinephase)
    
    %Creating the 2d sine wave
    %xp = x*cos(sinephase(si)) + y*sin(sinephase(si));
    %img = sin(2*pi*sinefreq*xp);
    
    %Creating the 2d sine wave
    img = [];
    des_len1 = time_points(1:length(time_points)*0.3); %Desired length for one kind of direction
    des_len2 = time_points(length(time_points)*0.3+1 : length(time_points));
    for ei = 1:num_electrodes
        temp_img1 = amplitude*sin(2*pi*sinefreq*des_len1 + (ei-1)*sinephase(si));
        temp_img2 = amplitude*sin(2*pi*sinefreq*des_len2 + -0.7*(ei-1)*sinephase(si));
        img(ei,:) = [temp_img1, temp_img2];
    end
    
    %Getting the 2d fft
    imgx = fftshift(fft2(img));
    powr2 = abs(imgx);
    %phas2 = angle(imgx);
    [FW, BW] = propagation_fun(img);

    %update plots
    set(imageh, 'CData', img);
    set(amph, 'CData', powr2);
    set(phaseh, 'YData', [FW, BW]);
    sgtitle('#2 kind of disruption')
    pause(1.0)
end
%} 

%% Function to calculate FW and BW propagation

function [thX1,thX2]=findingXboundaries(twod_fft,fx)

    %this is just a makiavelic way to find the boundaries of the interesting frequencies (i.e. thX1,thX2) ##to improve
    sumFFT=sum(twod_fft);
    threshold=0.05*max(sumFFT);
    [aa,var]=sort(abs(sumFFT-threshold));
    var=var(1:40); %40 big enough not to cut. ##to improve
    thX1=fx(min(var));
    thX2=fx(max(var));

end

function [FW, BW] = propagation_fun(egg_sig)
    
     twod_fft = abs(fftshift(fft2(egg_sig)));
     a=1; %width pixel to compute the frequency

     [m,n]=size(egg_sig); %numberPixel x and y

     fx=(1/a)*((1:n)-mean(1:n));
     fy=(1/a)*((1:m)-mean(1:m));

     [thX1,thX2]=findingXboundaries(twod_fft,fx);

     upperPart=twod_fft(fy>0,(fx>0 & fx<thX2));
     lowerPart=twod_fft(fy<0,(fx>0 & fx<thX2));

      FW = nanmean(upperPart(:));
      BW = nanmean(lowerPart(:));
end


%%
