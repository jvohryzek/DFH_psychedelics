%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Global DMT
% %
% % Jakub Vohryzek January 2026
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data info
clear all

subjects = {'01','02','03','06','07','09','10','11','12','13','15','17','18','19','22','23','25'} % all subject

numSbj = size(subjects,2);
parcelName = ['schaefer1000']; % schaefer1000
cnd = {'pcb_pre','pcb_post','dmt_pre','dmt_post'};
numCnd = size(cnd,2);
Tau=1;
NLAG=6;

%% FILTER SETTINGS
TR = 2;                 
fnq = 1/(2 * TR);               % Nyquist frequency
flp = 0.008; % 0.01;             % lowpass frequency of filter
fhi = 0.08; % 0.1;              % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter

%% Part 1. DATA
BOLD_processed_all.pcb_pre     = []; BOLD_processed_all.pcb_post     = [];
BOLD_processed_all.dmt_pre     = []; BOLD_processed_all.dmt_post     = [];

BOLD_processed_all_idx.pcb_pre = []; BOLD_processed_all_idx.pcb_post = [];
BOLD_processed_all_idx.dmt_pre = []; BOLD_processed_all_idx.dmt_post = [];

TS = cell(15,4); % for the optimisation part

% extract the timeseries for each subject
for c = 1:numCnd
    for subj = 1:size(subjects,2)
        prefix=['/Users/jakub/Datasets/Psychedelics/DMT/func_processed/','/sub-',subjects{subj},'/ses-',cnd{c},'/'];
        suffix=['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest_',parcelName,'.mat'];

        name = subjects{subj};
        disp(['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest'])
        file_mat = [prefix '/' suffix];
        XBOLD = cell(2,2,size(subjects,2)); % for the optimisation part
        % pre = 0; post = 0;

        BOLD = load(file_mat);
        TS{subj,c} = BOLD.tps;
        if strcmp(parcelName,'aal')
            order = [1:2:90 90:-2:1];
            [numAreas, numTp] = size(BOLD.tps(order,:)); % taking only aal90
        else
            [numAreas, numTp] = size(BOLD.tps);
        end
        Isubdiag = find(tril(ones(numAreas),-1));

        % BOLD signal to BOLD phase using the Hilbert transform
        excTp    = 0;                      % already taking 10 timepoints on line  87 ; number of timepoints to exclude from before and after
        BOLD_sbj{subj,c} = BOLD.tps;

        if strcmp(parcelName,'aal')
            order = [1:2:90 90:-2:1];
            [BOLD_processed, Phase_BOLD, staticFC] = ts2filtered(BOLD.tps(order,:), numAreas, bfilt, afilt, excTp);
        else
            [BOLD_processed, Phase_BOLD, staticFC] = ts2filtered(BOLD.tps, numAreas, bfilt, afilt, excTp);
        end
        % insideout reversibility
        ts = BOLD_processed(:,10:end-10);
        Tm = size(ts,2);
        FCtf = corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');       %% Core...FC tau foward
        FCtr = corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)'); %% FC tau reversal
        Itauf = -0.5*log(1-FCtf.*FCtf);  %% Mutual information
        Itaur = -0.5*log(1-FCtr.*FCtr);

        Itau_diff = abs(Itauf - Itaur);
        % Itau_diff(isnan(Itau_diff)) = eps; % accounting for missing regions
        Itau_diff_strength = sum(Itau_diff);
        Itau_diff_strength_sorted = sort(Itau_diff_strength,'descend');

        Reference = ((Itauf(Isubdiag)-Itaur(Isubdiag)).^2)';
        Reference2 = Itau_diff(Isubdiag);
        
        index = find(Reference>quantile(Reference,0)); % here 0, below I do it for the whole sweep
        
        FowRev = nanmean(Reference(index));
        FowRev_hierarchy = nanstd(Reference(index));
        FowRev_Fano = (FowRev_hierarchy.^2)/FowRev;

        %% calculating only a quantile of the difference of non-reversibility
        q_p_idx=0;
        for q_p = 0:0.05:0.95
            q_p_idx=q_p_idx+1;
            index_partial = find(Reference > quantile(Reference,q_p));
            FowRev_partial(q_p_idx)  = nanmean(Reference(index_partial));
            FowRev_std_partial(q_p_idx)  = std(Reference(index_partial));
            FowRev_Fano_partial(q_p_idx) = (FowRev_std_partial(q_p_idx).^2)/FowRev_partial(q_p_idx);
        end


    %% SAVING INTO DATA STRUCTURE

    % 1) BOLD for kmeans
    BOLD_processed_sbj.(cnd{c})(subj,:,:)           = BOLD_processed;
    BOLD_processed_all.(cnd{c})                     = cat(2,BOLD_processed_all.(cnd{c}),BOLD_processed);
    staticFC_sbj.(cnd{c})(subj,:,:)                 = staticFC;

    % idx
    BOLD_processed_sbj_idx.(cnd{c})(subj,:)         = str2num(subjects{subj})*ones(1,numTp);
    BOLD_processed_all_idx.(cnd{c})                 = cat(2,BOLD_processed_all_idx.(cnd{c}) ,str2num(subjects{subj})*ones(1,numTp));


    insideout_Fano.(cnd{c})(subj)                   = FowRev_Fano;

    end
end

%% figure non-reversibility
qnt = 20;
h1 = figure


% figure non-reversibility Fano

tmp_pcb_pre = insideout_Fano.pcb_pre;
tmp_pcb_post = insideout_Fano.pcb_post;
tmp_dmt_pre = insideout_Fano.dmt_pre;
tmp_dmt_post = insideout_Fano.dmt_post;

vp = violinplot([tmp_pcb_pre;tmp_pcb_post;tmp_dmt_pre;tmp_dmt_post]',{'PCB pre','PCB post','DMT pre','DMT post'});
vp(1).ViolinColor = [0.5 0.8 0.5];vp(2).ViolinColor = [0.5 0.8 0.5];vp(3).ViolinColor = [0.5 0.8 0.5];vp(4).ViolinColor = [0.5 0.5 0.8];

title(['DMT Fano - ',parcelName]); ylabel('Non-Reversibility Hierarchy')
% statistics
[hh1(1,1) p1(1,1)] = ttest(tmp_pcb_pre,tmp_pcb_post)
[hh1(1,2) p1(1,2)] = ttest(tmp_dmt_pre,tmp_dmt_post)
[hh1(1,3) p1(1,3)] = ttest(tmp_pcb_post,tmp_dmt_post)

if p1(1,1)<(0.05/(size(hh1,2)))
    plot(2,max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])+.001,'*r');hold on;
    plot([1 2],[max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004 max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004],'k','LineWidth',2)
elseif p1(1,1)<(0.05)
    plot(1.5,max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])+.001,'*k');hold on;
    plot([1 2],[max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004 max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004],'k','LineWidth',2)
end 
if p1(1,2)<(0.05/(size(hh1,2)))
    plot(3.5,max([mean(tmp_dmt_pre,2); mean(tmp_dmt_post,2)])+.001,'*r');hold on;
    plot([3 4],[max([mean(tmp_dmt_pre,2); mean(tmp_dmt_post,2)])-0.0008 max([mean(tmp_dmt_pre,2); mean(tmp_dmt_post,2)])-0.0008],'k','LineWidth',2)
elseif p1(1,2)<(0.05)
    plot(3.5,max([mean(tmp_dmt_pre,2); mean(tmp_dmt_post,2)])+.001,'*k');hold on;
    plot([3 4],[max([mean(tmp_dmt_pre,2); mean(tmp_dmt_post,2)])-0.0008 max([mean(tmp_dmt_pre,2); mean(tmp_dmt_post,2)])-0.0008],'k','LineWidth',2)
end 
if p1(1,3)<(0.05/(size(hh1,2)))
    plot(3,max([mean(tmp_pcb_post,2); mean(tmp_dmt_post,2)])+.001,'*r');hold on;
    plot([2 4],[max([mean(tmp_pcb_post,2); mean(tmp_dmt_post,2)])-0.0002 max([mean(tmp_pcb_post,2); mean(tmp_dmt_post,2)])-0.0002],'k','LineWidth',2)
elseif p1(1,3)<(0.05)
    plot(3,max([mean(tmp_pcb_post,2); mean(tmp_dmt_post,2)])+.001,'*k');hold on;
    plot([2 4],[max([mean(tmp_pcb_post,2); mean(tmp_dmt_post,2)])-0.0002 max([mean(tmp_pcb_post,2); mean(tmp_dmt_post,2)])-0.0002],'k','LineWidth',2)
end 


h1.Position = [100 100 1540 1400]
