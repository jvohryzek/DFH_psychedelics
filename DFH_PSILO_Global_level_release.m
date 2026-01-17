%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Global PSILO
% %
% % Jakub Vohryzek January 2026
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data info
clear all

subjects={'01','03','05','07','09','11','13','14','15'}
numSbj = size(subjects,2);
parcelName = 'schaefer1000'; 
cnd = {'pcb_rest1','pcb_rest2','psilo_rest1','psilo_rest2'};
numCnd = size(cnd,2);
Tau=1;

%% FILTER SETTINGS
TR = 3;                 
fnq = 1/(2 * TR);               % Nyquist frequency
flp = 0.008; % 0.01;             % lowpass frequency of filter
fhi = 0.08; % 0.1;              % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter
excTp    = 0;                      % number of timepoints to exclude from before and after


%% Part 1. DATA
BOLD_processed_all.pcb_rest1     = []; BOLD_processed_all.pcb_rest2     = [];
BOLD_processed_all.psilo_rest1     = []; BOLD_processed_all.psilo_rest2     = [];

BOLD_processed_all_idx.pcb_rest1 = []; BOLD_processed_all_idx.pcb_rest2 = [];
BOLD_processed_all_idx.psilo_rest1 = []; BOLD_processed_all_idx.psilo_rest2 = [];

% extract the timeseries for each subject
for c = 1:numCnd
    for subj = 1:size(subjects,2)
        prefix=['/Users/jakub/Datasets/Psychedelics/PSILO/func_processed/','/sub-',subjects{subj},'/ses-',cnd{c},'/'];
        suffix=['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest_',parcelName,'.mat'];

        name = subjects{subj};
        disp(['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest'])
        file_mat = [prefix '/' suffix];
        XBOLD = cell(2,2,size(subjects,2)); % for the optimisation part
        % pre = 0; post = 0;

        BOLD = load(file_mat);
        % XBOLD{2,2,size(subjects,2)} = BOLD;
        if strcmp(parcelName,'aal')
            [numAreas, numTp] = size(BOLD.tps(1:90,:)); % taking only aal90
        else
            [numAreas, numTp] = size(BOLD.tps); % taking only aal90

        end    
        Isubdiag = find(tril(ones(numAreas),-1));
        BOLD_sbj{subj,c} = BOLD.tps;

        % BOLD signal to BOLD phase using the Hilbert transform
        [BOLD_processed, Phase_BOLD, staticFC] = ts2filtered(BOLD.tps, numAreas, bfilt, afilt, excTp);

        % insideout reversibility
        ts = BOLD_processed(:,10:end-10);
        Tm = size(ts,2);
        FCtf = corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');       %% Core...FC tau foward
        FCtr = corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)'); %% FC tau reversal
        Itauf = -0.5*log(1-FCtf.*FCtf);  %% Mutual information... what is the this step for? Is it to work in positive values?
        Itaur = -0.5*log(1-FCtr.*FCtr);
        Itau_diff = abs(Itauf - Itaur);
        %Itau_diff(isnan(Itau_diff)) = eps; % accounting for missing regions
        Itau_diff_strength = sum(Itau_diff);
        Itau_diff_strength_sorted = sort(Itau_diff_strength,'descend');
        Reference = ((Itauf(Isubdiag)-Itaur(Isubdiag)).^2)';
        Reference2 = Itau_diff(Isubdiag);
        index = find(Reference>quantile(Reference,0)); % here 0, below I do it for the whole sweep
        FowRev = nanmean(Reference(index));
        FowRev_hierarchy = std(Reference(index));
        FowRev_Fano = (FowRev_hierarchy.^2)/FowRev;

    %% SAVING INTO DATA STRUCTURE

    % 1) BOLD for kmeans
    BOLD_processed_sbj.(cnd{c})(subj,:,:)           = BOLD_processed;
    BOLD_processed_all.(cnd{c})                     = cat(2,BOLD_processed_all.(cnd{c}),BOLD_processed);
    staticFC_sbj.(cnd{c})(subj,:,:)                 = staticFC;

    % idx
    BOLD_processed_sbj_idx.(cnd{c})(subj,:)         = str2num(subjects{subj})*ones(1,numTp);
    BOLD_processed_all_idx.(cnd{c})                 = cat(2,BOLD_processed_all_idx.(cnd{c}) ,str2num(subjects{subj})*ones(1,numTp));

    Phase_BOLD_all.(cnd{c})(subj,:,:)               = Phase_BOLD;

    insideout_Fano.(cnd{c})(subj)                   = FowRev_Fano;
    end
end

%% figure non-reversibility
h1 = figure

tmp_pcb_pre = insideout_Fano.pcb_rest1;
tmp_pcb_post = insideout_Fano.pcb_rest2;
tmp_psilo_pre = insideout_Fano.psilo_rest1;
tmp_psilo_post = insideout_Fano.psilo_rest2;

vp = violinplot([tmp_pcb_pre;tmp_pcb_post;tmp_psilo_pre;tmp_psilo_post]',{'PCB pre','PCB post','PSILO pre','PSILO post'});
vp(1).ViolinColor = [0.5 0.8 0.5];vp(2).ViolinColor = [0.5 0.8 0.5];vp(3).ViolinColor = [0.5 0.8 0.5];vp(4).ViolinColor = [0.5 0.5 0.8];

title(['PSILO Fano - ',parcelName]); ylabel('Non-Reversibility Hierarchy')
% statistics
[hh1(1,1) p1(1,1)] = ttest(tmp_pcb_pre,tmp_pcb_post)
[hh1(1,2) p1(1,2)] = ttest(tmp_psilo_pre,tmp_psilo_post)
[hh1(1,3) p1(1,3)] = ttest(tmp_pcb_post,tmp_psilo_post)

if p1(1,1)<(0.05/(size(hh1,2)))
    plot(2,max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])+.001,'*r');hold on;
    plot([1 2],[max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004 max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004],'k','LineWidth',2)
elseif p1(1,1)<(0.05)
    plot(1.5,max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])+.001,'*k');hold on;
    plot([1 2],[max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004 max([mean(tmp_pcb_pre,2); mean(tmp_pcb_post,2)])-0.0004],'k','LineWidth',2)
end 
if p1(1,2)<(0.05/(size(hh1,2)))
    plot(3.5,max([mean(tmp_psilo_pre,2); mean(tmp_psilo_post,2)])+.001,'*r');hold on;
    plot([3 4],[max([mean(tmp_psilo_pre,2); mean(tmp_psilo_post,2)])-0.0008 max([mean(tmp_psilo_pre,2); mean(tmp_psilo_post,2)])-0.0008],'k','LineWidth',2)
elseif p1(1,2)<(0.05)
    plot(3.5,max([mean(tmp_psilo_pre,2); mean(tmp_psilo_post,2)])+.001,'*k');hold on;
    plot([3 4],[max([mean(tmp_psilo_pre,2); mean(tmp_psilo_post,2)])-0.0008 max([mean(tmp_psilo_pre,2); mean(tmp_psilo_post,2)])-0.0008],'k','LineWidth',2)
end 
if p1(1,3)<(0.05/(size(hh1,2)))
    plot(3,max([mean(tmp_pcb_post,2); mean(tmp_psilo_post,2)])+.001,'*r');hold on;
    plot([2 4],[max([mean(tmp_pcb_post,2); mean(tmp_psilo_post,2)])-0.0002 max([mean(tmp_pcb_post,2); mean(tmp_psilo_post,2)])-0.0002],'k','LineWidth',2)
elseif p1(1,3)<(0.05)
    plot(3,max([mean(tmp_pcb_post,2); mean(tmp_psilo_post,2)])+.001,'*k');hold on;
    plot([2 4],[max([mean(tmp_pcb_post,2); mean(tmp_psilo_post,2)])-0.0002 max([mean(tmp_pcb_post,2); mean(tmp_psilo_post,2)])-0.0002],'k','LineWidth',2)
end 


h1.Position = [100 100 1540 1400]
