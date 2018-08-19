function A2_eegStimuliFeaturesRegressor

%% ================ load ordered stimuli labels for each subject ================
clear
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);

if ~exist('stimuli_labels_all_subj.mat')
    subj={'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9'...
        'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24'};
    % no subj 15, 16, >25% rejected trials
    for i=1:length(subj) %i=2
        
        load(subj{i});
        label{:,i}=label_new; % combine them into one file
        
    end
    
    stimuli_label=label;
    
    save stimuli_labels_all_subj stimuli_label;
else
    disp('...file exists, load the file')
    load stimuli_labels_all_subj
end
%% ================ load ordered EEG data for each subject ================

if ~exist('erp_data_all_subj.mat')
    softwarefolder1='C:/Users/CHUANJI/Documents/eeglab13_6_5b/';
    softwarefolder2='C:/Users/CHUANJI/Documents/fieldtrip-20170317/';
    parentfolder='F:/UPDATE/P_1_MIA/3_MIA_ERP/4_Result/ERP';
    fsep='/';
    chidrenfolder='CNT';
    
    addpath(softwarefolder1)
    addpath(softwarefolder2)
    cd(parentfolder);
    
    subj={'1_zhangjie' '2_caimengtong' '3_zhuhongtao' '4_caomengyuan' '5_yuezhen'...
        '6_wangpengcheng' '7_wangxiaoqing' '8_huangwenqiang' '9_suzijia' '10_sunguoyong'...
        '11_zhangweinan' '12_liuxinhe' '13_zhaoyiming' '14_qiuyanhong' '17_fengsijia'...
        '18_lanxi' '19_sundan' '20_ruanheng' '21_qimengchen' '22_liuchaoran' '23_sunyixuan' '24_guoxiaoqing'};
    
    for i=1:size(subj,2)
        
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        % load data
        EEG = pop_loadset('filename',[subj{i} '_synctrej.set'],'filepath',[parentfolder fsep subj{i} fsep chidrenfolder fsep]);
        EEG = eeg_checkset( EEG );
        temp{i} = EEG.data;
        
    end
    erp_data=temp;
    addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
    cd(addr)
    save erp_data_all_subj erp_data;
else
    disp('...file exists, load the file')
    load erp_data_all_subj
end
%% ================ load stimuli labels and features ================

addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);
disp('...load features for audiovisual trials')
load('features_av')
% 1 - hue; 2 - saturation; 3 - value; 4 - motion1; 5 - motion2; 6 - motion3; 7 - motion4; 8 - motion5; 9 - motion6; 10 - motion7; 11 - pitch; 12 - tempo; 13 - mode.
disp('...load features labels for audiovisual trials')

labels_avFile = [addr '\features_labels_av.txt'];
fileID = fopen(labels_avFile)
m=[];
m = textscan(fileID ,'%s','headerlines',0)
fclose(fileID)
disp('...labels for audiovisual trials were loaded')
labels_avOrder = m{1,1};

%---------------------------------------------
disp('...load features for video trials')
load('features_v')

labels_vFile = [addr '\features_labels_v.txt'];
fileID = fopen(labels_vFile)
m=[];
m = textscan(fileID ,'%s','headerlines',0)
fclose(fileID)
disp('...labels for video trials were loaded')
labels_vOrder = m{1,1};

%----------------------------------------------
disp('...load features for music trials')
load('features_a')

labels_aFile = [addr '\features_labels_a.txt'];
fileID = fopen(labels_aFile)
m=[];
m = textscan(fileID ,'%s','headerlines',0)
fclose(fileID)
disp('...labels for music trials were loaded')
labels_aOrder = m{1,1};
%% ================ get the stimuli labels, erp data just for accepted trials ================
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);
if ~exist('trials_reject.mat')
    
    softwarefolder1='C:/Users/CHUANJI/Documents/eeglab13_6_5b/';
    softwarefolder2='C:/Users/CHUANJI/Documents/fieldtrip-20170317/';
    parentfolder='F:/UPDATE/P_1_MIA/3_MIA_ERP/4_Result/ERP';
    fsep='/';
    chidrenfolder='CNT';

    addpath(softwarefolder1)
    addpath(softwarefolder2)
    cd(parentfolder);

    subj={'1_zhangjie' '2_caimengtong' '3_zhuhongtao' '4_caomengyuan' '5_yuezhen'...
        '6_wangpengcheng' '7_wangxiaoqing' '8_huangwenqiang' '9_suzijia' '10_sunguoyong'...
        '11_zhangweinan' '12_liuxinhe' '13_zhaoyiming' '14_qiuyanhong' '17_fengsijia'...
        '18_lanxi' '19_sundan' '20_ruanheng' '21_qimengchen' '22_liuchaoran' '23_sunyixuan' '24_guoxiaoqing'};

    for i=1:size(subj,2)

        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        % load data
        EEG = pop_loadset('filename',[subj{i} '_synctrej.set'],'filepath',[parentfolder fsep subj{i} fsep chidrenfolder fsep]);
        EEG = eeg_checkset( EEG );
        trials_reject{i} = EEG.reject.rejmanual;

    end

    addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
    cd(addr);
    save trials_reject trials_reject
else
    disp('...file exists, load the file')
    load trials_reject
end

if ~exist('erp_data_accepted.mat')
    addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
    cd(addr);
    
    subj={'1_zhangjie' '2_caimengtong' '3_zhuhongtao' '4_caomengyuan' '5_yuezhen'...
        '6_wangpengcheng' '7_wangxiaoqing' '8_huangwenqiang' '9_suzijia' '10_sunguoyong'...
        '11_zhangweinan' '12_liuxinhe' '13_zhaoyiming' '14_qiuyanhong' '17_fengsijia'...
        '18_lanxi' '19_sundan' '20_ruanheng' '21_qimengchen' '22_liuchaoran' '23_sunyixuan' '24_guoxiaoqing'};
    
    for ni=1:size(subj,2)
        ii=1;
        for i=1:length(stimuli_label{1,ni})
            if trials_reject{1,ni}(i)==0;
                stimuli_label_accepted{ii,ni}=stimuli_label{1,ni}(i,:);
                erp_data_accepted{1,ni}(:,:,ii)=erp_data{1,ni}(:,:,i);
                ii=ii+1;
            end
        end
    end   
    save stimuli_label_accepted stimuli_label_accepted
    save erp_data_accepted erp_data_accepted
else
    disp('...file exists, load the file')
    load stimuli_label_accepted
    load erp_data_accepted
end

if ~exist('erp_data_av.mat')
    for ni=1:length(erp_data_accepted) % subject
        ii=1;
        for i=1:size(erp_data_accepted{1,ni},3) % trials
            if length(stimuli_label_accepted{i,ni}{1,1})==7 % select audiovisual trials
                stimuli_label_av{ii,ni}=stimuli_label_accepted{i,ni}{1,1};
                erp_data_av{1,ni}(:,:,ii)=erp_data_accepted{1,ni}(:,:,i);
                ii=ii+1;
            end
        end
    end
    save erp_data_av erp_data_av
    save stimuli_label_av stimuli_label_av
else
    disp('...file exists, load the file')
    load stimuli_label_av
    load erp_data_av
end

if ~exist('erp_data_v.mat')
    for ni=1:length(erp_data_accepted) % subject
        ii=1;
        for i=1:size(erp_data_accepted{1,ni},3) % trials
            if (length(stimuli_label_accepted{i,ni}{1,1})==5) && (strcmp(stimuli_label_accepted{i,ni}{1,1}(1),'V')) % select audiovisual trials
                stimuli_label_v{ii,ni}=stimuli_label_accepted{i,ni}{1,1};
                erp_data_v{1,ni}(:,:,ii)=erp_data_accepted{1,ni}(:,:,i);
                ii=ii+1;
            end
        end
    end
    save erp_data_v erp_data_v
    save stimuli_label_v stimuli_label_v
else
    disp('...file exists, load the file')
    load stimuli_label_v
    load erp_data_v
end

if ~exist('erp_data_a.mat')
    for ni=1:length(erp_data_accepted) % subject
        ii=1;
        for i=1:size(erp_data_accepted{1,ni},3) % trials
            if (length(stimuli_label_accepted{i,ni}{1,1})==5) && (strcmp(stimuli_label_accepted{i,ni}{1,1}(1),'M')) % select audiovisual trials
                stimuli_label_a{ii,ni}=stimuli_label_accepted{i,ni}{1,1};
                erp_data_a{1,ni}(:,:,ii)=erp_data_accepted{1,ni}(:,:,i);
                ii=ii+1;
            end
        end
    end
    save erp_data_a erp_data_a
    save stimuli_label_a stimuli_label_a
else
    disp('...file exists, load the file')
    load stimuli_label_a
    load erp_data_a
end
%% ================== match the features with erp data =====================
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);

if ~exist('features_av_ordered.mat')
    
    for ni = 1:size(stimuli_label_av,2)
        i=1;
        parfor i=1:size(erp_data_av{1,ni},3) % trials
            ind=[];
            ind = find(strcmp(stimuli_label_av{i,1},labels_avOrder));
            disp('...getting the ordered features of audiovisual trials')
            features_av_ordered(i,ni) = features_av(ind,:);
        end
    end
    load features_av_ordered features_av_ordered
else
    disp('...file exists, load the file')
    load features_av_ordered  
end




