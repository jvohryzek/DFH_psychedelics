clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % LZ Complexity LSD
% %
% % Jakub Vohryzek January 2026
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data info


subjects={'01','02','03','04','06','09','10','11','12','13','15','17','18','19','20'}

numSbj = size(subjects,2);
parcelName = {'schaefer100','schaefer200','schaefer300','schaefer400','schaefer500','schaefer1000'}; % schaefer1000
cnd = {'pcb_rest1','pcb_rest2','pcb_rest3','lsd_rest1','lsd_rest2','lsd_rest3'};
numCnd = size(cnd,2);
Tau=1;
%% FILTER SETTINGS
TR = 2;                 
fnq = 1/(2 * TR);               % Nyquist frequency
flp = 0.04; % 0.01;             % lowpass frequency of filter
fhi = 0.07; % 0.1;              % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter

% extract the timeseries for each subject
for parc = 1:size(parcelName,2)
    %% Part 1. DATA
    BOLD_processed_all.pcb_rest1     = []; BOLD_processed_all.pcb_rest2     = []; BOLD_processed_all.pcb_rest3     = [];
    BOLD_processed_all.lsd_rest1     = []; BOLD_processed_all.lsd_rest2     = []; BOLD_processed_all.lsd_rest3     = [];

    BOLD_processed_all_idx.pcb_rest1 = []; BOLD_processed_all_idx.pcb_rest2 = []; BOLD_processed_all_idx.pcb_rest3     = [];
    BOLD_processed_all_idx.lsd_rest1 = []; BOLD_processed_all_idx.lsd_rest2 = []; BOLD_processed_all_idx.lsd_rest3     = [];


    
    TS = cell(15,6); % for the optimisation part

    for c = 1:numCnd
        for subj = 1:size(subjects,2)

            prefix=['/Users/jakub/Datasets/Psychedelics/LSD/func_processed/','/sub-',subjects{subj},'/ses-',cnd{c},'/'];
            suffix=['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest_',parcelName{parc},'.mat'];
    
            name = subjects{subj};
            disp(['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest'])
            file_mat = [prefix '/' suffix];
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
            Itauf = -0.5*log(1-FCtf.*FCtf);  %% Mutual information... what is the this step for? Is it to work in positive values?
            Itaur = -0.5*log(1-FCtr.*FCtr);
    
            Itau_diff = abs(Itauf - Itaur);
            Itau_diff_strength = sum(Itau_diff);
            Itau_diff_strength_sorted = sort(Itau_diff_strength,'descend');
            x_temp = [0:0.01:1];

            Reference = ((Itauf(Isubdiag)-Itaur(Isubdiag)).^2)';
            Reference2 = Itau_diff(Isubdiag);
            
            index = find(Reference>quantile(Reference,0)); % here 0, below I do it for the whole sweep
            
            FowRev = nanmean(Reference(index));
            FowRev_hierarchy = std(Reference(index));
            FowRev_Fano = (FowRev_hierarchy.^2)/FowRev;
            FowRev_CV = FowRev_hierarchy/FowRev;

            insideout.(cnd{c})(subj,parc)                        = FowRev;
            insideout_hierarchy.(cnd{c})(subj,parc)              = FowRev_hierarchy;
            insideout_Fano.(cnd{c})(subj,parc)                   = FowRev_Fano;
            insideout_CV.(cnd{c})(subj,parc)                     = FowRev_CV;
            index_hierarchy.(cnd{c}){subj,parc}                  = index;

            
            [C_norm ,C] = LZcomplexity(zscore(BOLD_processed,0,2), numAreas);
            LZcom.(cnd{c})(subj,parc)                            = mean(C_norm);
            %LZcom_median.(cnd{c})(subj,parc)                     = median(C_norm);
            clear C_nom C



        end
    end
end
%%
insideout_Fano_pcb = (insideout_Fano.pcb_rest1+insideout_Fano.pcb_rest3)./2;
insideout_Fano_lsd = (insideout_Fano.lsd_rest1+insideout_Fano.lsd_rest3)./2;

LZcom_pcb = (LZcom.pcb_rest1+LZcom.pcb_rest3)./2;
LZcom_lsd = (LZcom.lsd_rest1+LZcom.lsd_rest3)./2;

insideout_Fano_pcb1 = (insideout_Fano.pcb_rest1)
insideout_Fano_lsd1 = (insideout_Fano.lsd_rest1)

LZcom_pcb1 = (LZcom.pcb_rest1);
LZcom_lsd1 = (LZcom.lsd_rest1);

insideout_Fano_pcb3 = (insideout_Fano.pcb_rest3);
insideout_Fano_lsd3 = (insideout_Fano.lsd_rest3);

LZcom_pcb3 = (LZcom.pcb_rest3);
LZcom_lsd3 = (LZcom.lsd_rest3);
%% Figure


h1 = figure
i=4
    tmp_lsd_pre  = LZcom_pcb(:,i)';
    tmp_lsd_post = LZcom_lsd(:,i)';

    vp = violinplot([tmp_lsd_pre;tmp_lsd_post]',{'PCB pre','PCB post','DMT pre','DMT post'});
    vp(1).ViolinColor = [0.5 0.8 0.5];vp(2).ViolinColor = [0.5 0.5 0.8];
    
    title(['LSD LZcom ',parcelName{i}]); ylabel('LSD LZcom')
    % statistics
    [hh1(1,i) p1(1,i)] = ttest(tmp_lsd_pre,tmp_lsd_post)

    
    if p1(1,i)<(0.05/(size(hh1,2)))
        plot(2,max([mean(tmp_lsd_pre,2); mean(tmp_lsd_post,2)])+.001,'*r');hold on;
        plot([1 2],[max([mean(tmp_lsd_pre,2); mean(tmp_lsd_post,2)])-0.0004 max([mean(tmp_lsd_pre,2); mean(tmp_lsd_post,2)])-0.0004],'k','LineWidth',2)
    elseif p1(1,i)<(0.05)
        plot(1.5,max([mean(tmp_lsd_pre,2); mean(tmp_lsd_post,2)])+.001,'*k');hold on;
        plot([1 2],[max([mean(tmp_lsd_pre,2); mean(tmp_lsd_post,2)])-0.0004 max([mean(tmp_lsd_pre,2); mean(tmp_lsd_post,2)])-0.0004],'k','LineWidth',2)
    end 


