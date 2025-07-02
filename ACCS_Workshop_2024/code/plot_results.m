%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The author of the code is Sanket S Houde working under Dr.Pragathi
% Balasubramani at Department of Cognitive Sciences, IIT Kanpur. You can
% reach me or Dr. Pragathi via our email IDs. 
% Sanket S Houde - sankethoude@gmail.com
% Dr. Pragathi Balasubramani - pbalasub@iitk.ac.in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('C:\Users\sanke\Downloads\Workshop_Sanket\PAC_values');
%addpath('C:\Users\sanke\Downloads\SANKET WORKSHOP\PAC_values');

%% Gastric PAC 
pac_gastric_theta = [];
pac_gastric_alpha = [];
pac_gastric_beta = [];

task1 = load('PAC_T1_new_workshop.mat');
task2 = load('PAC_T2_new_workshop.mat');

pac_gastric_theta(:,1) = task1.pac_gastric_top5(1,:);
pac_gastric_theta(:,2) = task2.pac_gastric_top5(1,:);

pac_gastric_alpha(:,1) = task1.pac_gastric_top5(2,:);
pac_gastric_alpha(:,2) = task2.pac_gastric_top5(2,:);

pac_gastric_beta(:,1) = task1.pac_gastric_top5(3,:);
pac_gastric_beta(:,2) = task2.pac_gastric_top5(3,:);

freq = {'Frontal', 'Occipital', 'Parietal','Central','Left temporal','Right temporal'};
figure;
bar(freq,pac_gastric_theta);
title('Theta gastric PAC');
ylim([0,0.3]);
legend({'Eyes Open Task','Breathe in Breathe out Task'});

figure;
bar(freq,pac_gastric_alpha);
title('Alpha gastric PAC');
ylim([0,0.3]);
legend({'Eyes Open Task','Breathe in Breathe out Task'});

figure;
bar(freq,pac_gastric_beta);
title('Beta gastric PAC');
ylim([0,0.3]);
legend({'Eyes Open Task','Breathe in Breathe out Task'});

%% Heart PAC
pac_heart_theta = [];
pac_heart_alpha = [];
pac_heart_beta = [];

pac_heart_theta(:,1) = task1.pac_heart_top5(1,:);
pac_heart_theta(:,2) = task2.pac_heart_top5(1,:);

pac_heart_alpha(:,1) = task1.pac_heart_top5(2,:);
pac_heart_alpha(:,2) = task2.pac_heart_top5(2,:);

pac_heart_beta(:,1) = task1.pac_heart_top5(3,:);
pac_heart_beta(:,2) = task2.pac_heart_top5(3,:);

freq = {'Frontal', 'Occipital', 'Parietal','Central','Left temporal','Right temporal'};
figure;
bar(freq,pac_heart_theta);
title('Theta Heart PAC');
ylim([0,0.14]);
legend({'Eyes Open Task','Breathe in Breathe out Task'});

figure;
bar(freq,pac_heart_alpha);
title('Alpha Heart PAC');
ylim([0,0.14]);
legend({'Eyes Open Task','Breathe in Breathe out Task'});

figure;
bar(freq,pac_heart_beta);
title('Beta Heart PAC');
ylim([0,0.14]);
legend({'Eyes Open Task','Breathe in Breathe out Task'});

