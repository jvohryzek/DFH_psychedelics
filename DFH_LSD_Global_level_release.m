%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Global LSD
% %
% % Jakub Vohryzek January 2026
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data info
clear all;close all;

subjects={'01','02','03','04','06','09','10','11','12','13','15','17','18','19','20'}
numSbj = size(subjects,2);
parcelName = ['schaefer1000']; 
cnd = {'pcb_rest1','pcb_rest2','pcb_rest3','lsd_rest1','lsd_rest2','lsd_rest3'};
numCnd = size(cnd,2);
Tau=1;

%% FILTER SETTINGS
TR = 2;                 
fnq = 1/(2 * TR);               % Nyquist frequency
flp = 0.008; % 0.01;             % lowpass frequency of filter
fhi = 0.08; % 0.1;              % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter
excTp    = 0;                      % number of timepoints to exclude from before and after

%% Part 1. DATA
BOLD_processed_all.pcb_rest1     = []; BOLD_processed_all.pcb_rest2     = []; BOLD_processed_all.pcb_rest3     = [];
BOLD_processed_all.lsd_rest1     = []; BOLD_processed_all.lsd_rest2     = []; BOLD_processed_all.lsd_rest3     = [];

BOLD_processed_all_idx.pcb_rest1 = []; BOLD_processed_all_idx.pcb_rest2 = []; BOLD_processed_all_idx.pcb_rest3     = [];
BOLD_processed_all_idx.lsd_rest1 = []; BOLD_processed_all_idx.lsd_rest2 = []; BOLD_processed_all_idx.lsd_rest3     = [];
XBOLD = cell(15,4); % for the optimisation part

% extract the timeseries for each subject
for c = 1:numCnd
    for subj = 1:size(subjects,2)
        prefix=['/Users/jakub/Datasets/Psychedelics/LSD/func_processed/','/sub-',subjects{subj},'/ses-',cnd{c},'/'];
        suffix=['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest_',parcelName,'.mat'];

        name = subjects{subj};
        disp(['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest'])
        file_mat = [prefix '/' suffix];
        % pre = 0; post = 0;

        BOLD = load(file_mat);
        if c == 1 
             XBOLD{subj,1} = BOLD.tps;
        elseif c == 3
             XBOLD{subj,2} = BOLD.tps;
        elseif c == 4
             XBOLD{subj,3} = BOLD.tps;
        elseif c == 6
             XBOLD{subj,4} = BOLD.tps;
        end
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
        Itau_diff(isnan(Itau_diff)) = eps; % accounting for missing regions
        Itau_diff_strength = sum(Itau_diff);
        Itau_diff_strength_sorted = sort(Itau_diff_strength,'descend');
        Reference = ((Itauf(Isubdiag)-Itaur(Isubdiag)).^2)';
        Reference2 = Itau_diff(Isubdiag);
        index = find(Reference>quantile(Reference,0)); % here 0, below I do it for the whole sweep
        FowRev = nanmean(Reference(index));
        FowRev_hierarchy = std(Reference(index));

        FowRev_Fano = (FowRev_hierarchy.^2)/FowRev;

        %% calculating only a quantile of the difference of non-reversibilitty
        q_p_idx=0;
        for q_p = 0:0.05:0.95
            q_p_idx=q_p_idx+1;
            index_partial = find(Reference > quantile(Reference,q_p));
            FowRev_partial(q_p_idx)  = nanmean(Reference(index_partial));
            FowRev_integrated = sum(FowRev_partial)./size(FowRev_partial,2);
            FowRev_std_partial(q_p_idx)  = std(Reference(index_partial));
            FowRev_Fano_partial(q_p_idx) = (FowRev_std_partial(q_p_idx).^2)/FowRev_partial(q_p_idx);
        end
        %% calculating gradients
        %[coe,cc,lat,tss,exp,mu] = pca(Itau_diff);


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

%%
figure


tmp_pcb_rest = (insideout_Fano.pcb_rest1+insideout_Fano.pcb_rest3)./2;
tmp_lsd_rest = (insideout_Fano.lsd_rest1+insideout_Fano.lsd_rest3)./2;

vp = violinplot([tmp_pcb_rest;tmp_lsd_rest]',{'PCB','LSD'});
vp(1).ViolinColor = [0.5 0.8 0.5];vp(2).ViolinColor = [0.5 0.5 0.8];

title(['LSD Fano Temporal Asymmetry ',parcelName]); ylabel('Fano Temporal Asymmetry')
% statistics
[hh3 p3] = ttest(tmp_pcb_rest,tmp_lsd_rest)



