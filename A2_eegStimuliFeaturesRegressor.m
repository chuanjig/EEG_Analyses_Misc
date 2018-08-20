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
features_av_ordered = NaN(268,13,22);
if ~exist('features_av_ordered.mat')
    
    for ni = 1:size(stimuli_label_av,2) % ni = 22
        i=1;
        for i=1:size(erp_data_av{1,ni},3) % trials
            ind=[];
            ind = find(strcmp(stimuli_label_av{i,ni},labels_avOrder));
            features_av_ordered(i,:,ni) = features_av(ind,:);
        end
        disp('...getting the ordered features of audiovisual trials')
        disp(['...participant ' num2str(ni) ' is done'])
    end
    save features_av_ordered features_av_ordered
else
    disp('...file exists, load the file')
    load features_av_ordered  
end

features_v_ordered = NaN(54,10,22);
if ~exist('features_v_ordered.mat')
    
    for ni = 1:size(stimuli_label_v,2) % ni = 22
        i=1;
        for i=1:size(erp_data_v{1,ni},3) % trials
            ind=[];
            ind = find(strcmp(stimuli_label_v{i,ni},labels_vOrder));
            features_v_ordered(i,:,ni) = features_v(ind,:);
        end
        disp('...getting the ordered features of audiovisual trials')
        disp(['...participant ' num2str(ni) ' is done'])
    end
    save features_v_ordered features_v_ordered
else
    disp('...file exists, load the file')
    load features_v_ordered  
end

features_a_ordered = NaN(53,3,22);
if ~exist('features_a_ordered.mat')
    
    for ni = 1:size(stimuli_label_a,2) % ni = 22
        i=1;
        for i=1:size(erp_data_a{1,ni},3) % trials
            ind=[];
            ind = find(strcmp(stimuli_label_a{i,ni},labels_aOrder));
            features_a_ordered(i,:,ni) = features_a(ind,:);
        end
        disp('...getting the ordered features of audiovisual trials')
        disp(['...participant ' num2str(ni) ' is done'])
    end
    save features_a_ordered features_a_ordered
else
    disp('...file exists, load the file')
    load features_a_ordered  
end

disp('...now deleting some electrodes that we do not need')

if ~exist('erp_data_av_new.mat')
    for ni = 1:size(erp_data_av,2) % number of participants
        disp('...doing it for the audiovisual trials')
        temp=erp_data_av{1,ni};
        temp([33 43 65 66 67 68],:,:)=[];
        erp_data_av_new{1,ni}=temp;
        clear temp
        disp('...doing it for the visual only trials')
        temp=erp_data_v{1,ni};
        temp([33 43 65 66 67 68],:,:)=[];
        erp_data_v_new{1,ni}=temp;
        clear temp
        disp('...doing it for the auditory only trials')
        temp=erp_data_a{1,ni};
        temp([33 43 65 66 67 68],:,:)=[];
        erp_data_a_new{1,ni}=temp;
        clear temp
    end
    save erp_data_av_new erp_data_av_new
    save erp_data_v_new erp_data_v_new
    save erp_data_a_new erp_data_a_new
else
    disp('...file exists, load the file')
    load erp_data_av_new
    load erp_data_v_new
    load erp_data_a_new
end

%% ================== doing regression =====================
disp('...now all of the order of the features are matched with ERP data')
disp('...now we could begin to do the regression on the ERP data');

disp('...regressing out saturation for unimodal video')
disp('...then check out the results')

% 1 - hue; 2 - saturation; 3 - value; 4 - motion1; 5 - motion2; 6 - motion3; 7 - motion4; 8 - motion5; 9 - motion6; 10 - motion7; 11 - pitch; 12 - tempo; 13 - mode.
% 1 - hue; 2 - saturation; 3 - value; 4 - motion1; 5 - motion2; 6 - motion3; 7 - motion4; 8 - motion5; 9 - motion6; 10 - motion7;
% 1 - pitch; 2 - tempo; 3 - mode.
if ~exist('erp_data_v_final.mat')
    ni = [];
    for ni = 1:size(erp_data_v_new,2) % ni = 2
        erp_data_v_ing = [];
        features_v_ing = [];
        
        erp_data_v_ing = erp_data_v_new{1,ni};
        erp_data_v_ing = reshape(erp_data_v_ing, [size(erp_data_v_ing,1)*size(erp_data_v_ing,2), size(erp_data_v_ing,3)]);
        erp_data_v_ing = erp_data_v_ing';
        disp('...choosing the feature to be regressed, saturation')
        features_v_ing = features_v_ordered(1:size(erp_data_v_ing,1),2,ni);
        
        disp('...normalize before regression')
        ifeatures_v_ing = [];
        regmat = [];
        for ifeatures_v_ing = 1:length(features_v_ing)
            regmat(ifeatures_v_ing) = log10(features_v_ing(ifeatures_v_ing));
        end
        regmat=regmat';
        resid = [];
        regmat(:,end+1)=1; % add the additional column here.
        disp('...doing regression')
        for j=1:size(erp_data_v_ing,2)
            [B1,BINT1,R1] = regress(erp_data_v_ing(:,j),regmat);
            resid(:,j)=R1;
        end
        resid=resid';
        disp('...assigning new values')
        erp_data_v_new{1,ni} = reshape(resid, [size(erp_data_v_new{1,ni},1) size(erp_data_v_new{1,ni},2) size(erp_data_v_new{1,ni},3)]);
        disp(['subj ' num2str(ni) ' is done'])
    end
    disp('...saving data')
    erp_data_v_final = erp_data_v_new;
    save erp_data_v_final erp_data_v_final
    disp('...regressing out saturation is done')
else
    disp('...file exists, load the file')
    load erp_data_v_final   
end



