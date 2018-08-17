function A2_eegStimuliFeaturesRegressor

%% ================ load ordered stimuli labels for each subject ================
clear
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);
subj={'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9'...
    'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24'};
% no subj 15, 16, >25% rejected trials
for i=1:length(subj) %i=2
    
    load(subj{i});
    label{:,i}=label_new; % combine them into one file
    
end

stimuli_label=label;

save stimuli_labels_all_subj stimuli_label;
load stimuli_labels_all_subj
%% ================ load ordered EEG data for each subject ================

clear

softwarefolder1='C:/Users/CHUANJI/Documents/eeglab13_6_5b/';
softwarefolder2='C:/Users/CHUANJI/Documents/fieldtrip-20170317/';
parentfolder='F:/UPDATE/P_1_MIA/3_MIA_ERP/4_Result/ERP';
fsep='/';
chidrenfolder='CNT'

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
cd(addr);
save erp_data_all_subj erp_data;
clear
load erp_data_all_subj

%% ================ run PCA for motion parameters ================
clear,clc
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures';
cd(addr);
load('features') 
% hue saturation value motion1 motion2 motion3 motion4 motion5 motion6 motion7 pitch tempo 
% video_valence music_valence video_arousal music_arousal condition_vector1 condition_vector2 condition_vector3
temp=features(:,4:10); % extract all 7 motion features
[coeff,score,latent,tsquared,explained,mu] = pca(temp);
% explained = 98.3953, the first component explain a lot, I would just used the first component.
motion_pca=temp*coeff; % multiply Eigenvectors(coeff) and raw data, you get the data, the first column is the first component.
motion_pca=motion_pca(:,1); % get the first component of motion features.

%% ================ update the feature matrix with pca-ed motion =================
features_pca=features(:,1:3);
features_pca(:,4)=motion_pca;
features_pca(:,5:6)=features(:,11:12); % hue saturation brightness motion pitch tempo (6 columns)
features_pca(:,7:13) =features(:,13:19); % video_valence music_valence video_arousal music_arousal condition_vector1 condition_vector2 condition_vector3
save features_pca features_pca

%% ================ get the labels and features from excel and combined them into one matrix manually ===============
clear,clc
addr='F:\UPDATE\P_1_MIA\3_MIA_ERP\5_Analysis\EEG\erp_forward_encoding_model';
cd(addr);
load stimuli_labels_formatted % 270*1
load features_pca % 270*6
% I just copy data from features_excel into label_excel and save it as stimuli_labels_features (270*7).

%% ================ get the stimuli labels, erp data just for accepted trials============

clear

softwarefolder1='/Users/Gao/Documents/MATLAB/eeglab13_6_5b/';
softwarefolder2='/Users/Gao/Documents/MATLAB/fieldtrip-master/';
parentfolder='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/4_Result/ERP';
fsep='/';
chidrenfolder='CNT'

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

addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';
cd(addr);
save trials_reject trials_reject

clear,clc

subj={'1_zhangjie' '2_caimengtong' '3_zhuhongtao' '4_caomengyuan' '5_yuezhen'...
    '6_wangpengcheng' '7_wangxiaoqing' '8_huangwenqiang' '9_suzijia' '10_sunguoyong'...
    '11_zhangweinan' '12_liuxinhe' '13_zhaoyiming' '14_qiuyanhong' '17_fengsijia'...
    '18_lanxi' '19_sundan' '20_ruanheng' '21_qimengchen' '22_liuchaoran' '23_sunyixuan' '24_guoxiaoqing'};

addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';
cd(addr);

load trials_reject trials_reject
load stimuli_labels_all_subj % stimuli labels for 450 trials all subjects (22)
load erp_data_all_subj % erp data for 450 trials all subjects (22)

for ni=1:size(subj,2)
    ii=1;
    for i=1:length(stimuli_label{1,ni})
        if trials_reject{1,ni}(i)==0;
            temp_label{ii,ni}=stimuli_label{1,ni}(i,:);
            temp_data{1,ni}(:,:,ii)=erp_data{1,ni}(:,:,i);
            ii=ii+1;
        end
    end 
end

% ================ get the stimuli labels, erp data just for audiovisual trials============

clear i ii ni trials_reject stimuli_label erp_data stimuli_label_new erp_data_new
for ni=1:length(temp_data) % subject
    ii=1;
    for i=1:size(temp_data{1,ni},3) % trials
        if length(temp_label{i,ni}{1,1})==7 % delete unimodal trials
            stimuli_label_new{ii,ni}=temp_label{i,ni}{1,1};
            erp_data_new{1,ni}(:,:,ii)=temp_data{1,ni}(:,:,i);
            ii=ii+1;
        end
    end
end

save stimuli_labels_all_subj_new stimuli_label_new
save erp_data_all_subj_new erp_data_new

%% ================== match the features with erp data =====================
clear,clc
addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';

cd(addr);
load stimuli_labels_all_subj_new  % labels
load erp_data_all_subj_new % erp data
load stimuli_labels_features % features
load features_pca

% clear features_pca_ordered
for ni = 1:size(erp_data_new,2) % number of participants
    
    for i=1:size(erp_data_new{1,ni},3) % number of trials
        for ii=1:270
            if stimuli_label_new{i,ni}==stimuli_labels_features{ii,1}
                
                features_pca_new(i,1)=stimuli_labels_features{ii,2}; % hue
                features_pca_new(i,2)=stimuli_labels_features{ii,3}; % saturation
                features_pca_new(i,3)=stimuli_labels_features{ii,4}; % value
                features_pca_new(i,4)=stimuli_labels_features{ii,5}; % motion
                features_pca_new(i,5)=stimuli_labels_features{ii,6}; % pitch
                features_pca_new(i,6)=stimuli_labels_features{ii,7}; % tempo
                features_pca_new(i,7)=stimuli_labels_features{ii,8}; 
                features_pca_new(i,8)=stimuli_labels_features{ii,9};
                features_pca_new(i,9)=stimuli_labels_features{ii,10};
                features_pca_new(i,10)=stimuli_labels_features{ii,11};
                features_pca_new(i,11)=stimuli_labels_features{ii,12};
                features_pca_new(i,12)=stimuli_labels_features{ii,13};
                features_pca_new(i,13)=stimuli_labels_features{ii,14};

            end
        end
    end
    
    features_pca_ordered{ni}=features_pca_new;
    disp(['Running ============= ** participant ' num2str(ni) ' ======== is FINISHED'])
    clear features_pca_new
end

save features_pca_ordered features_pca_ordered

%% ================== get ready all the data needed =====================
clear,clc
addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';

cd(addr);

load features_pca_ordered % ordered features
load erp_data_all_subj_new % ordered erp data

% ---------- get z-transformed attribute values
for ni = 1:size(erp_data_new,2) % number of participants
    
    for j=1:size(features_pca_ordered{1,ni},2)
        temp=features_pca_ordered{1,ni}(:,j);
        stdvalue=std(temp);
        meanvalue=mean(temp);
        for i=1:size(temp,1)
            features_pca_ordered_new(i,j)= (temp(i)- meanvalue)./stdvalue;
        end
    end
    features_pca_ordered_z_trans{ni}=features_pca_ordered_new;
    disp(['Running ============= ** participant ' num2str(ni) ' ======== is FINISHED'])
    clear features_pca_ordered_new temp
end

save features_pca_ordered_z_trans features_pca_ordered_z_trans

% ---------- delete electrodes: 33-M1, 43-M2, 65-HEO, 66-VEO, 67-GFP, 68-REF.

load erp_data_all_subj_new % final ordered erp data

for ni = 1:size(erp_data_new,2) % number of participants
    
    temp=erp_data_new{1,ni};
    temp([33 43 65 66 67 68],:,:)=[];
    erp_data_new_elec{1,ni}=temp;
    clear temp
    
end

save erp_data_all_subj_new_elec erp_data_new_elec

% ---------- get the 500-900ms time window (350points to 550points) data.
% ---------- get the 150-250ms time window (175points to 225points) data.

clear,clc
load erp_data_all_subj_new_elec

for ni = 1:size(erp_data_new_elec,2) % number of participants
    
    temp=erp_data_new_elec{1,ni};
%     temp1=temp(:,[351:550],:);
    temp1=temp(:,[175:225],:);
    temp2=mean(temp1, 2); % average across time points
    temp3=squeeze(temp2);
    temp4=permute(temp3, [2 1]); % permute the data to trial by electrode
    erp_data_new_elec_n200{ni}=temp4;
    clear temp temp1 temp2 temp3 temp4
    
end
save erp_data_all_subj_new_elec_n200 erp_data_new_elec_n200
% save erp_data_all_subj_new_elec_lpp erp_data_new_elec_lpp

%% ================== step1: generate the AMs =====================

% ---------- do regression
clear,clc
load features_pca_ordered_z_trans
load erp_data_all_subj_new_elec_n200
% load erp_data_all_subj_new_elec_lpp

for ni=1:length(erp_data_new_elec_n200) % participant ni=2

    multiregmat=features_pca_ordered_z_trans{1, ni};
    multiregmat=multiregmat(:,1:6);
    multiregmat(:,end+1)=1; % add the additional column here.
    
    temp=erp_data_new_elec_n200{1, ni};
    
    nTrain = round(size(temp,1)./6)*5; % The 5/6 part was taken as the training data.
    
    %     size(temp(1:nTrain,:),1)
    %     size(temp(nTrain+1:size(temp,1),:),1)
    
    for j=1:62 % electrode
        [B1,BINT1,R1] = regress(temp(1:nTrain,j), multiregmat(1:nTrain,:));
        beta_coef(:,j,ni)=B1;
    end
    clear multiregmat temp B1 BINT1 R1

end
    
beta_coef=beta_coef(1:6,:,:);
beta_coef=permute(beta_coef, [2 1 3]); % electrode, 6 features, participants

save beta_lower_features beta_coef

% ---------- make attribute maps - using fieldtrip here.
clear,clc
softwarefolder='/Users/Gao/Documents/MATLAB/fieldtrip-master/';
addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';
addpath(softwarefolder)
cd(addr);
ft_defaults;

load eeg2fieldtrip_data_format
load beta_lower_features
data.label([33 43 65 66 67 68])=[]; % delete useless elect
data.elec.pnt([33 43 65 66 67 68],:)=[];
data.elec.label=data.label';
data.var([33 43 65 66 67 68],:)=[];
beta_coef=mean(beta_coef,3); % average across 22 participants
name={'Hue' 'Saturation' 'Brightness' 'Motion' 'Pitch' 'Tempo'};
fontnamevalue='Arial';
fontsizevalue=15;

for i=1:size(beta_coef, 2)
    
    subplot(2,3,i);
    data.avg=repmat(beta_coef(:,i), 1, 600); % make repeated 600 columns, 62*600
    cfg = [];
    cfg.parameter = 'avg';
    cfg.ylim = [-1.5 1.5];
    cfg.zlim = [-1.5 1.5];
    cfg.comment   = 'no';
    cfg.style     = 'both';
    % cfg.style     = 'fill';
    cfg.layout    = 'quickcap64.mat';
    ft_topoplotER(cfg, data)
    title(name{i}, 'FontName', fontnamevalue, 'FontSize', fontsizevalue)
    
    %     address='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/6_Figures/foward_encoding_model';
    %     cd(address)
    %     set(gcf, 'PaperPositionMode', 'auto');
    %     print([name{i} '_attribute'],'-dpng','-r300');
    %     close;
end

c=colorbar('Ticks',[-1.5,0,1.5], 'FontSize', 10);
set(c,'position',[.05 .4 .03 .3])

address='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/6_Figures/foward_encoding_model';
cd(address)
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 attribute_maps.eps
print('attribute_maps','-dpng','-r300');
close;

%% ================== step2: generate the PMs =====================
% in the 1/6 part = features * attribute maps, get observed maps

clear,clc
addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';
cd(addr);
load beta_lower_features
load features_pca_ordered_z_trans % features
load erp_data_all_subj_new_elec_n200
% load erp_data_all_subj_new_elec_lpp

% beta_coef=mean(beta_coef, 3); % average across 22 participants
beta_coef=permute(beta_coef, [2 1 3]); % 6*62

for ni=1:size(features_pca_ordered_z_trans,2)
    beta_temp=[];
    beta_temp=beta_coef(:,:,ni);
    temp=features_pca_ordered_z_trans{1,ni};
    nTrain = round(size(temp,1)./6)*5; % The 5/6 part was taken as the training data.
    
    predicted_map{ni} = temp(nTrain+1:size(temp,1),1:6)*beta_temp; % get predicted maps

    observed_map{ni} = erp_data_new_elec_n200{1,ni}(nTrain+1:size(temp,1),:); % get observed maps
    
end

%% ================== step3: running correlation between AMs and PMs =====================
% ---------- run correlation between predicted maps and observed maps for each stimulus
for ni=1:length(observed_map)
    for i=1:size(observed_map{1,ni},1)
        for j=1:size(predicted_map{1,ni},1)
            [r,p] = corrcoef(observed_map{1,ni}(i,:),predicted_map{1,ni}(j,:),'rows','pairwise');
            crR(i,j) = r(1,2);
            sigP(i,j) = p(1,2);
        end
    end
    crR_all{ni}=crR;
    sigP_all{ni}=sigP;
    clear crR sigP
end
save corr_sig crR_all sigP_all
% crR, sigP, are nTest*nTest matrices.
% rows: actual (observed)
% columns: predicted

% ---------- plot correlation map based on crR

clims=[-1 1];
fig=imagesc(crR_all{1,2}, clims);

set(gca, 'ydir', 'normal');
axis square;

set(gcf,'Color','w')

% add the colorbar, make it prettier
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Label.String = 'pearson R';
drawnow;

axpos = get(gca, 'Position');
set(gca,'position',[axpos(1),axpos(2),axpos(3),axpos(4)]);
drawnow;

hTitle  = title ('Correlation Strength (Pearson R)', 'FontName', 'Arial', 'FontSize', 20);
hXLabel = xlabel('Predicted Labels' ,'FontName', 'Arial', 'FontSize', 15);
hYLabel = ylabel('Observed Labels','FontName', 'Arial', 'FontSize', 15);

set(hTitle, 'units', 'normalized');
tiPos = get (hTitle, 'position'); % axis x, y value of title position
xyrange = axis;

set(hTitle, 'position', [0.52 1.02 0]); % move title up a bit

set(gca, ...
    'FontName'   , 'Arial'    , ...
    'color'      , 'none'     , ...
    'xgrid'      , 'off'      ,...
    'ygrid'      , 'off'      ,...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'FontSize'   ,  15, ...
    'LineWidth'   , 1         );

address='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/6_Figures/foward_encoding_model';
cd(address)
set(gcf, 'PaperPositionMode', 'auto');
print('correlation_map','-dpng','-r300');
close;

% -------------- find the rank of the crR of the predicted (column) for each observed (row)
clear,clc
addr='/Volumes/Gao/UPDATE/P_1_MIA/3_MIA_ERP/5_Analysis/EEG/erp_forward_encoding_model';
cd(addr);
load corr_sig

% clear prediction_acc
for ni=1:22
temp=crR_all{1,ni};
for i=1:size(temp,1)
    [B Ind]=sort(temp(i,:),'descend') ;
    index=find(Ind==i);
    prediction_acc(i)=(size(temp,1)-index)./(size(temp,1)-1);
end
prediction_all{ni}=prediction_acc;
clear temp prediction_acc
end

for ni=1:22
    
    accuracy(ni)=mean(prediction_all{ni})
    
end
    
mean(accuracy)
max(accuracy)
min(accuracy)

