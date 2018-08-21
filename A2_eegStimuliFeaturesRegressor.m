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

disp('==================================================================')
disp('...regressing out saturation for unimodal video')
disp('...then check out the results')
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);

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

disp('==================================================================')
disp('...regressing out saturation for unimodal music')
disp('...then check out the results')
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);

% 1 - hue; 2 - saturation; 3 - value; 4 - motion1; 5 - motion2; 6 - motion3; 7 - motion4; 8 - motion5; 9 - motion6; 10 - motion7; 11 - pitch; 12 - tempo; 13 - mode.
% 1 - hue; 2 - saturation; 3 - value; 4 - motion1; 5 - motion2; 6 - motion3; 7 - motion4; 8 - motion5; 9 - motion6; 10 - motion7;
% 1 - pitch; 2 - tempo; 3 - mode.
if ~exist('erp_data_a_final.mat')
    ni = [];
    for ni = 1:size(erp_data_a_new,2) % ni = 2
        erp_data_a_ing = [];
        features_a_ing = [];
        
        erp_data_a_ing = erp_data_a_new{1,ni};
        erp_data_a_ing = reshape(erp_data_a_ing, [size(erp_data_a_ing,1)*size(erp_data_a_ing,2), size(erp_data_a_ing,3)]);
        erp_data_a_ing = erp_data_a_ing';
        disp('...choosing the feature to be regressed, pitch and mode')
        features_a_ing = features_a_ordered(1:size(erp_data_a_ing,1),[1 3],ni);
        
        disp('...normalize before regression')
        ifeatures_a_ing = [];
        regmat = [];
        for ifeatures_a_ing = 1:length(features_a_ing)
            regmat(ifeatures_a_ing) = log10(features_a_ing(ifeatures_a_ing));
        end
        regmat=regmat';
        resid = [];
        regmat(:,end+1)=1; % add the additional column here.
        disp('...doing regression')
        for j=1:size(erp_data_a_ing,2)
            [B1,BINT1,R1] = regress(erp_data_a_ing(:,j),regmat);
            resid(:,j)=R1;
        end
        resid=resid';
        disp('...assigning new values')
        erp_data_a_new{1,ni} = reshape(resid, [size(erp_data_a_new{1,ni},1) size(erp_data_a_new{1,ni},2) size(erp_data_a_new{1,ni},3)]);
        disp(['subj ' num2str(ni) ' is done'])
    end
    disp('...saving data')
    erp_data_a_final = erp_data_a_new;
    save erp_data_a_final erp_data_a_final
    disp('...regressing out pitch and mode is done')
else
    disp('...file exists, load the file')
    load erp_data_a_final
end


%% ================== get ready for ANOVA =====================
disp('==================================================================')
disp('...get ready for ANOVA for unimodal visual stimuli')

addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);
load erp_data_v_final

disp('...get indices for each condition')
ni = [];
VPTrial=[];
VXTrial=[];
VNTrial=[];
for ni = 1:size(erp_data_v_final,2) % ni = 2
    erp_data_v_final_ing = [];
    erp_data_v_final_ing = erp_data_v_final{1,ni};
    iTrial = [];
    for iTrial = 1:size(erp_data_v_final{1,ni},3) % iTrial = 2
        if strcmp(stimuli_label_v{iTrial,ni}(1:2),'VP')
            VPTrial(iTrial,ni) = iTrial;
        elseif strcmp(stimuli_label_v{iTrial,ni}(1:2),'VX')
            VXTrial(iTrial,ni) = iTrial;
        elseif strcmp(stimuli_label_v{iTrial,ni}(1:2),'VN')
            VNTrial(iTrial,ni) = iTrial;
        else
            disp('...not a trial')
        end
    end
end

disp('now got the indices for each condition')
disp('now start to get the mean amplitudes')
erp_data_VP=[];
erp_data_VX=[];
erp_data_VN=[];
for ni = 1:size(erp_data_v_final,2) % ni = 2
    erp_data_v_final_ing = [];
    erp_data_v_final_ing = erp_data_v_final{1,ni};
    disp('deleting 0s in a index vector')
    VPTrialTemp = [];
    VPTrialTemp = VPTrial(:,ni);
    VPTrialTemp(VPTrialTemp==0) = [];
    
    VXTrialTemp = [];
    VXTrialTemp = VXTrial(:,ni);
    VXTrialTemp(VXTrialTemp==0) = [];

    VNTrialTemp = [];
    VNTrialTemp = VNTrial(:,ni);
    VNTrialTemp(VNTrialTemp==0) = [];
    
    disp('get the mean amplitudes for condition')
    erp_data_VP{1,ni} = mean(erp_data_v_final_ing(:,:,VPTrialTemp),3);
    erp_data_VX{1,ni} = mean(erp_data_v_final_ing(:,:,VXTrialTemp),3);
    erp_data_VN{1,ni} = mean(erp_data_v_final_ing(:,:,VNTrialTemp),3);
end

% anterior: F5,F3,F1,FZ,F2,F4,F6; FC5,FC3,FC1,FCZ,FC2,FC4,FC6
% posterior: C5,C3,C1,CZ,C2,C4,C6; CP5,CP3,CP1,CPZ,CP2,CP4,CP6;P5,P3,P1,PZ,P2,P4,P6;PO5,PO3,POZ,PO4,PO6
% ---------- get the 100 to 150ms time window (150points to 175points) data. N125.
% ---------- get the 150 to 200ms time window (175points to 200points) data. P170.
% ---------- get the 220 to 320ms time window (210points to 260points) data. N250.
% ---------- get the 500 to 900ms time window (350points to 550points) data. LPP.
TimeWindow1 = 150:175; 
TimeWindow2 = 175:200; 
TimeWindow3 = 210:260; 
TimeWindow4 = 350:550;
anteriorInd = [7:13 16:22];
posteriorInd = [25:31 34:40 43:49 52:56];

disp('...processing N125')
N125_ANT_VP = [];
N125_POS_VP = [];
N125_ANT_VX = [];
N125_POS_VX = [];
N125_ANT_VN = [];
N125_POS_VN = [];

for ni = 1:size(erp_data_v_final,2) % ni = 2
    erp_data_VP_temp = [];
    erp_data_VP_temp = erp_data_VP{1,ni};
    erp_data_VX_temp = [];
    erp_data_VX_temp = erp_data_VX{1,ni};    
    erp_data_VN_temp = [];
    erp_data_VN_temp = erp_data_VN{1,ni};
    
    N125_ANT_VP(ni) = mean2(erp_data_VP_temp(anteriorInd,TimeWindow1));
    N125_POS_VP(ni) = mean2(erp_data_VP_temp(posteriorInd,TimeWindow1));
    N125_ANT_VX(ni) = mean2(erp_data_VX_temp(anteriorInd,TimeWindow1));
    N125_POS_VX(ni) = mean2(erp_data_VX_temp(posteriorInd,TimeWindow1));
    N125_ANT_VN(ni) = mean2(erp_data_VN_temp(anteriorInd,TimeWindow1));
    N125_POS_VN(ni) = mean2(erp_data_VN_temp(posteriorInd,TimeWindow1));    
  
end

disp('...processing P170')
P170_ANT_VP = [];
P170_POS_VP = [];
P170_ANT_VX = [];
P170_POS_VX = [];
P170_ANT_VN = [];
P170_POS_VN = [];

for ni = 1:size(erp_data_v_final,2) % ni = 2
    erp_data_VP_temp = [];
    erp_data_VP_temp = erp_data_VP{1,ni};
    erp_data_VX_temp = [];
    erp_data_VX_temp = erp_data_VX{1,ni};    
    erp_data_VN_temp = [];
    erp_data_VN_temp = erp_data_VN{1,ni};
    
    P170_ANT_VP(ni) = mean2(erp_data_VP_temp(anteriorInd,TimeWindow2));
    P170_POS_VP(ni) = mean2(erp_data_VP_temp(posteriorInd,TimeWindow2));
    P170_ANT_VX(ni) = mean2(erp_data_VX_temp(anteriorInd,TimeWindow2));
    P170_POS_VX(ni) = mean2(erp_data_VX_temp(posteriorInd,TimeWindow2));
    P170_ANT_VN(ni) = mean2(erp_data_VN_temp(anteriorInd,TimeWindow2));
    P170_POS_VN(ni) = mean2(erp_data_VN_temp(posteriorInd,TimeWindow2));    
  
end

disp('...processing N250')
N250_ANT_VP = [];
N250_POS_VP = [];
N250_ANT_VX = [];
N250_POS_VX = [];
N250_ANT_VN = [];
N250_POS_VN = [];

for ni = 1:size(erp_data_v_final,2) % ni = 2
    erp_data_VP_temp = [];
    erp_data_VP_temp = erp_data_VP{1,ni};
    erp_data_VX_temp = [];
    erp_data_VX_temp = erp_data_VX{1,ni};    
    erp_data_VN_temp = [];
    erp_data_VN_temp = erp_data_VN{1,ni};
    
    N250_ANT_VP(ni) = mean2(erp_data_VP_temp(anteriorInd,TimeWindow3));
    N250_POS_VP(ni) = mean2(erp_data_VP_temp(posteriorInd,TimeWindow3));
    N250_ANT_VX(ni) = mean2(erp_data_VX_temp(anteriorInd,TimeWindow3));
    N250_POS_VX(ni) = mean2(erp_data_VX_temp(posteriorInd,TimeWindow3));
    N250_ANT_VN(ni) = mean2(erp_data_VN_temp(anteriorInd,TimeWindow3));
    N250_POS_VN(ni) = mean2(erp_data_VN_temp(posteriorInd,TimeWindow3));    
  
end

disp('...processing LPP')
LPP_ANT_VP = [];
LPP_POS_VP = [];
LPP_ANT_VX = [];
LPP_POS_VX = [];
LPP_ANT_VN = [];
LPP_POS_VN = [];

for ni = 1:size(erp_data_v_final,2) % ni = 2
    erp_data_VP_temp = [];
    erp_data_VP_temp = erp_data_VP{1,ni};
    erp_data_VX_temp = [];
    erp_data_VX_temp = erp_data_VX{1,ni};    
    erp_data_VN_temp = [];
    erp_data_VN_temp = erp_data_VN{1,ni};
    
    LPP_ANT_VP(ni) = mean2(erp_data_VP_temp(anteriorInd,TimeWindow4));
    LPP_POS_VP(ni) = mean2(erp_data_VP_temp(posteriorInd,TimeWindow4));
    LPP_ANT_VX(ni) = mean2(erp_data_VX_temp(anteriorInd,TimeWindow4));
    LPP_POS_VX(ni) = mean2(erp_data_VX_temp(posteriorInd,TimeWindow4));
    LPP_ANT_VN(ni) = mean2(erp_data_VN_temp(anteriorInd,TimeWindow4));
    LPP_POS_VN(ni) = mean2(erp_data_VN_temp(posteriorInd,TimeWindow4));    
  
end

disp('...get the data matrix of VP ready')
N125 = [N125_ANT_VP' N125_POS_VP' N125_ANT_VX' N125_POS_VX' N125_ANT_VN' N125_POS_VN'];
P170 = [P170_ANT_VP' P170_POS_VP' P170_ANT_VX' P170_POS_VX' P170_ANT_VN' P170_POS_VN'];
N250 = [N250_ANT_VP' N250_POS_VP' N250_ANT_VX' N250_POS_VX' N250_ANT_VN' N250_POS_VN'];
LPP = [LPP_ANT_VP' LPP_POS_VP' LPP_ANT_VX' LPP_POS_VX' LPP_ANT_VN' LPP_POS_VN'];

unimodalV_RegressionFinal = [N125 P170 N250 LPP];
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);
save unimodalV_RegressionFinal unimodalV_RegressionFinal

