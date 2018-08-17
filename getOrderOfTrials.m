function getOrderOfTrials
% modified based on bahavior_analysis_model to deal with the ordre of trials used
% for regressing out the lower level features.

clear,clc
dir='F:\UPDATE\P_1_MIA\3_MIA_ERP\4_Result\BEHAVIOR';
cd (dir)

subj={'1_zhangjie' '2_caimengtong' '3_zhuhongtao' '4_caomengyuan' '5_yuezhen'...
    '6_wangpengcheng' '7_wangxiaoqing' '8_huangwenqiang' '9_suzijia' '10_sunguoyong'...
    '11_zhangweinan' '12_liuxinhe' '13_zhaoyiming' '14_qiuyanhong' '15_cidanzhuoga' '16_panzijian' '17_fengsijia'...
    '18_lanxi' '19_sundan' '20_ruanheng' '21_qimengchen' '22_liuchaoran' '23_sunyixuan' '24_guoxiaoqing'};

for i=1:8 %length(subj) %i=16

    clear m1 m2 m3 m4 m5 m6 %num1 num2 num3 label1 label2 label3 pre1 pre2 pre3 pre_v pre_a label_new label1_new label2_new label3_new pre1_v pre2_v pre3_v pre1_a pre2_a pre3_a % this is important otherwise dimension does not match.
mm1=[]; 
mm2=[];
mm3=[];
mm4=[];
mm5=[];
mm6=[];

%read all data from file into 'm' as strings with space delimiter
[m1(:,1),m1(:,2),m1(:,3),m1(:,4),m1(:,5),m1(:,6), m1(:,7), m1(:,8), m1(:,9),...
m1(:,10), m1(:,11), m1(:,12)] = textread([subj{i} '/exp1.log'],'%s %s %s %s %s %s %s %s %s %s %s %s')
m1(1:4,:)=[];
mm1=str2double(m1(:,3)); % 3rd column

[m2(:,1),m2(:,2),m2(:,3),m2(:,4),m2(:,5),m2(:,6), m2(:,7), m2(:,8), m2(:,9),...
m2(:,10), m2(:,11), m2(:,12)] = textread([subj{i} '/exp2.log'],'%s %s %s %s %s %s %s %s %s %s %s %s')
m2(1:4,:)=[];
mm2=str2double(m2(:,3)); % 3rd column

[m3(:,1),m3(:,2),m3(:,3),m3(:,4),m3(:,5),m3(:,6), m3(:,7), m3(:,8), m3(:,9),...
m3(:,10), m3(:,11), m3(:,12)] = textread([subj{i} '/exp3.log'],'%s %s %s %s %s %s %s %s %s %s %s %s')
m3(1:4,:)=[];
mm3=str2double(m3(:,3)); % 3rd column

[m4(:,1),m4(:,2),m4(:,3),m4(:,4),m4(:,5),m4(:,6), m4(:,7), m4(:,8), m4(:,9),...
m4(:,10), m4(:,11), m4(:,12)] = textread([subj{i} '/exp4.log'],'%s %s %s %s %s %s %s %s %s %s %s %s')
m4(1:4,:)=[];
mm4=str2double(m4(:,3)); % 3rd column

[m5(:,1),m5(:,2),m5(:,3),m5(:,4),m5(:,5),m5(:,6), m5(:,7), m5(:,8), m5(:,9),...
m5(:,10), m5(:,11), m5(:,12)] = textread([subj{i} '/exp5.log'],'%s %s %s %s %s %s %s %s %s %s %s %s')
m5(1:4,:)=[];
mm5=str2double(m5(:,3)); % 3rd column

[m6(:,1),m6(:,2),m6(:,3),m6(:,4),m6(:,5),m6(:,6), m6(:,7), m6(:,8), m6(:,9),...
m6(:,10), m6(:,11), m6(:,12)] = textread([subj{i} '/exp6.log'],'%s %s %s %s %s %s %s %s %s %s %s %s')
m6(1:4,:)=[];
mm6=str2double(m6(:,3)); % 3rd column

% ----subject
currDir = 'F:\UPDATE\P_1_MIA\3_MIA_ERP\8_Paper\1_ERPs_paper\Submission\BiologicalPsychology\ResponseStimuliFeatures\';
clear label1_new label2_new label3_new label1 label2 label3 pre1_v pre2_v pre3_v pre1_a pre2_a pre3_a
filename = [currDir 'SUBJECT' num2str(i) '.xlsx'];
[num1 label1] = xlsread(filename,1);
[num2 label2] = xlsread(filename,2);
[num3 label3] = xlsread(filename,3);

pre1=[mm1;mm2]; % unique code from 1 to 150;
pre2=[mm3;mm4];
pre3=[mm5;mm6];

% pre1 ----
for ni=1:length(pre1)
    if pre1(ni) < 200
        pre1_v(ni)=pre1(ni+1)-200; %valence
        pre1_a(ni)=pre1(ni+2)-200; %arousal
    end
end
        ind_pre1_v=find(pre1_v~=0);
        pre1_v=pre1_v(ind_pre1_v)'; % valence values 150 * 1
        ind_pre1_a=find(pre1_a~=0);
        pre1_a=pre1_a(ind_pre1_a)'; % arousal values 150 * 1
        ind_pre1=find(pre1<200);
        pre1_new=pre1(ind_pre1); % code (in the order of presentation) 150*1
        
for mt=1:length(pre1_new)
    label1_new{mt}=label1{pre1_new(mt)};
end
label1_new=label1_new';  % stimuli labels (in the order of presentation) 150*1

ni=[];

% pre2 ----
for ni=1:length(pre2)
    if pre2(ni) < 200
        pre2_v(ni)=pre2(ni+1)-200; %valence
        pre2_a(ni)=pre2(ni+2)-200; %arousal
    end
end
        ind_pre2_v=find(pre2_v~=0);
        pre2_v=pre2_v(ind_pre2_v)'; % valence values 150 * 1
        ind_pre2_a=find(pre2_a~=0);
        pre2_a=pre2_a(ind_pre2_a)'; % arousal values 150 * 1
        ind_pre2=find(pre2<200);
        pre2_new=pre2(ind_pre2); % code (in order of presentation) 150*1
mt=[];  
for mt=1:length(pre2_new)
    label2_new{mt}=label2{pre2_new(mt)};
end
label2_new=label2_new'; % stimuli labels (in the order of presentation) 150*1

% pre3 ----
ni=[];
for ni=1:length(pre3)
    if pre3(ni) < 200
        pre3_v(ni)=pre3(ni+1)-200; %valence
        pre3_a(ni)=pre3(ni+2)-200; %arousal
    end
end
        ind_pre3_v=find(pre3_v~=0);
        pre3_v=pre3_v(ind_pre3_v)'; % valence values 150 * 1
        ind_pre3_a=find(pre3_a~=0);
        pre3_a=pre3_a(ind_pre3_a)'; % arousal values 150 * 1
        ind_pre3=find(pre3<200);
        pre3_new=pre3(ind_pre3); % code (in order of presentation) 150*1
mt=[];  
for mt=1:length(pre3_new)
    label3_new{mt}=label3{pre3_new(mt)};
end
label3_new=label3_new'; % stimuli labels (in the order of presentation) 150*1

label_new=[label1_new;label2_new;label3_new]; %stimuli labels for all trials in one subject
pre_v=[pre1_v;pre2_v;pre3_v]; %valence ratings
pre_a=[pre1_a;pre2_a;pre3_a]; %arousal ratings

% save
save([currDir 'subj' num2str(i) '.mat'], 'label_new', 'pre_v', 'pre_a')
end

