%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Graining Global PSILO
% %
% % Jakub Vohryzek January 2026
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data info
clear all;


subjects={'01','03','05','07','09','11','13','14','15'}

numSbj = size(subjects,2);
parcelName = {'schaefer100','schaefer200','schaefer300','schaefer400','schaefer500','schaefer1000'}; % schaefer1000
cnd = {'pcb_rest1','pcb_rest2','psilo_rest1','psilo_rest2'};
numCnd = size(cnd,2);
Tau=1;
NLAG=6;

%% FILTER SETTINGS
TR = 3;                 
fnq = 1/(2 * TR);               % Nyquist frequency
flp = 0.008; % 0.01;             % lowpass frequency of filter
fhi = 0.08; % 0.1;              % highpass
Wn = [flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k = 2;                          % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn);   % construct the filter

% extract the timeseries for each subject
for parc = 1:size(parcelName,2)
    %% Part 1. DATA
    BOLD_processed_all.pcb_rest1     = []; BOLD_processed_all.pcb_rest2     = [];
    BOLD_processed_all.psilo_rest1     = []; BOLD_processed_all.psilo_rest2     = [];
    
    BOLD_processed_all_idx.pcb_rest1 = []; BOLD_processed_all_idx.pcb_rest2 = [];
    BOLD_processed_all_idx.psilo_rest1 = []; BOLD_processed_all_idx.psilo_rest2 = [];
        
    TS = cell(9,4); % for the optimisation part

    for c = 1:numCnd
        for subj = 1:size(subjects,2)
        prefix=['/Users/jakub/Datasets/Psychedelics/PSILO/func_processed/','/sub-',subjects{subj},'/ses-',cnd{c},'/'];
        suffix=['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest_',parcelName{parc},'.mat'];

        name = subjects{subj};
        disp(['sub-',subjects{subj},'_ses-',cnd{c},'_task-rest'])
        file_mat = [prefix '/' suffix];
        XBOLD = cell(2,2,size(subjects,2)); % for the optimisation part
        % pre = 0; post = 0;
    
            BOLD = load(file_mat);
            TS{subj,c} = BOLD.tps;
            if strcmp(parcelName{parc},'aal')
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

        insideout_Fano.(cnd{c})(subj,parc)                   = FowRev_Fano;

        end
    end
end

%% statistics
figure
subplot(2,1,1)
violinplot(abs(insideout_Fano.(cnd{4})-insideout_Fano.(cnd{3})))
ylabel('Post-Pre Temporal Asymmetry');xlabel('Parcelation Resolution');
ylim([0,0.025])

for i = 1:6
    [h_DMT(i) p_DMT(i)] = ttest(insideout_Fano.(cnd{3})(:,i),insideout_Fano.(cnd{4})(:,i))
    [h_PCB(i) p_PCB(i)] = ttest(insideout_Fano.(cnd{1})(:,i),insideout_Fano.(cnd{2})(:,i))

end
subplot(2,1,2)
plot(p_DMT);ylabel('P-values');xlabel('Parcelation Resolution');%ylim([0, 0.45])
hold on
plot([1 6],[0.05 0.05],'r--');xlim([1,6])


