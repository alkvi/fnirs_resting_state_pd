% Subject-level correlation analysis
% Ran on R2024a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/resting_state

%% Load BIDS

% Which session are we computing for?
session_to_run = "right";

% Load the NIRx probe used
probe_folder = "data/nirx_probe_" + session_to_run;
nirx_probe = nirs.io.loadNIRxProbe(probe_folder,true);

% Load BIDS
my_data_dir = 'data/park_move_rs_fnirs_bids';
raw_data = nirs.bids.loadBIDS(my_data_dir, true, nirx_probe);

% Where to store data
matrix_folder = "data/subject_matrices";

% Bad short channels
bad_shorts_data = readtable("data/bad_shorts.csv");

% See what we loaded
demographics = nirs.createDemographicsTable(raw_data);
disp(demographics);

%% Save raw data

%save('data/mat_files/raw_data_rsfnirs.mat','raw_data');

% or load
%raw_data = importdata('data/mat_files/raw_data.mat');

%% Visualize probe

raw_data_test = raw_data(1);
raw_data_test.probe.defaultdrawfcn='?';

figure();
raw_data_test.probe.defaultdrawfcn='3D label mes h';
raw_data_test.probe.draw;

%% Add age to demographics

age_data = readtable('data/basic_demographics.csv'); 
demographics.age = NaN(height(demographics),1);

% Add age data
for idx=1:height(age_data)
    subj_id_seek = string(age_data.subject(idx)); 
    match_idx = strcmp(replace(demographics.bids_subject, "rs", ""), subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).age = ones(sum(match_idx),1) * age_data(idx,:).age;
end

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
raw_data = job.run(raw_data);

%% Mark bad short channels in raw_data

% We will go through each subject and mark the 
% probe table with bad channels
for subj_idx=1:length(raw_data)
    
    % Find matching subject/session for data lookup
    subj_id_seek = string(raw_data(subj_idx).demographics.bids_subject); 
    session_seek = string(raw_data(subj_idx).demographics.session);
    match_idx = strcmp(bad_shorts_data.subject, subj_id_seek) ...
        & strcmp(bad_shorts_data.session, session_seek);
    if sum(match_idx) < 1
        disp("ERROR: no match")
        pause;
    end
    
    % Add a column 'bad' in probe.link
    raw_data(subj_idx).probe.link.bad(:) = 0;

    % Get bad shorts for this subject
    bad_shorts = string(bad_shorts_data(match_idx,:).bad_shorts{1});
    if (bad_shorts == "")
        continue
    end
    
    % Format is e.g. S8-D11
    % go through and mark bad entries
    bad_list = split(bad_shorts,',');
    for bad_list_idx=1:length(bad_list)
        src_det = bad_list(bad_list_idx);
        src_det = split(src_det, "-");
        src = str2num(erase(src_det(1), "S"));
        det = str2num(erase(src_det(2), "D"));
        link = raw_data(subj_idx).probe.link;
        link_idx = (link.source == src) & (link.detector == det);
        
        raw_data(subj_idx).probe.link.bad(link_idx) = 1;
        fprintf("Marked src-det %d-%d for subj %s session %s as bad\n", src, det, subj_id_seek, session_seek);
    end

end

%% Run on all subjects

for subj_idx=1:length(raw_data)

    % Some subject strings for saving data
    subject_id = raw_data(subj_idx).demographics.bids_subject;
    session = raw_data(subj_idx).demographics.session;

    % Just run one side for now
    if ~strcmp(session, session_to_run)
        continue;
    end

    % Trim the data
    d = raw_data(subj_idx).data;
    t = raw_data(subj_idx).time;
    valid_time = t >= 10 & t <= max(t) - 10;
    d_trim = d(valid_time,:);
    t_trim = t(valid_time,:);
    raw_data(subj_idx).data = d_trim;
    raw_data(subj_idx).time = t_trim;

    % Run individual analysis
    fprintf("Getting ConnStats for for subj %s session %s\n", subject_id, session);
    [ConnStats, Graph, hb] = get_connectivity_stats(raw_data(subj_idx));
    
    % The actual matrices from the analysis
    pearson_matrix = ConnStats.R;
    adjacency_matrix = Graph.adjacency;

    % Remove entries in matrices for short channels
    short_separation_idx = (ConnStats.probe.link.detector >= 16);
    short_separation_idx = find(short_separation_idx > 0);
    short_separation_idx = sort(short_separation_idx,'descend');
    link = ConnStats.probe.link;
    for i=1:length(short_separation_idx)
        ss_idx = short_separation_idx(i);
        pearson_matrix(ss_idx,:) = [];
        pearson_matrix(:,ss_idx) = [];
        adjacency_matrix(ss_idx,:) = [];
        adjacency_matrix(:,ss_idx) = [];
        link(ss_idx,:) = [];
    end

    % Also remove bad quality channels
    % S4-D3, S4-D4, around S8
    bad_src = [4, 4, 8, 8, 8, 8];
    bad_det = [3, 4, 6, 7, 8, 9];
    for pair_idx=1:length(bad_src)
        s = bad_src(pair_idx);
        d = bad_det(pair_idx);
        bad_idx = (link.source == s) & (link.detector == d);

        pearson_matrix(bad_idx,:) = [];
        pearson_matrix(:,bad_idx) = [];
        adjacency_matrix(bad_idx,:) = [];
        adjacency_matrix(:,bad_idx) = [];
        link(bad_idx,:) = [];
    end
    
    % Write the final matrices to file
    filebase = matrix_folder + "/" + string(subj_idx) + "_" + subject_id + "_" + session;
    writematrix(pearson_matrix, filebase + "_R.xlsx")
    writematrix(adjacency_matrix, filebase + "_adjacency.xlsx")
end

%% Get connectivity

function [ConnStats, Graph, hb] = get_connectivity_stats(raw_data)

% Label which channels are short-separation.
job = nirs.modules.LabelShortSeperation();
job.max_distance = 8;
raw_data = job.run(raw_data);

% Some sources will connect to adjacent short-distance detectors which will
% look like a long channel but actually is not. Mark these also.
for i=1:size(raw_data,1)
    short_separation_idx =(raw_data(i).probe.link.detector >= 16);
    raw_data(i).probe.link.ShortSeperation(short_separation_idx) = 1;
end

% Run MBLL
j=nirs.modules.OpticalDensity;
j=nirs.modules.BeerLambertLaw(j);
hb = j.run(raw_data);

% Only look at HbO.
job=nirs.modules.KeepTypes;
job.types={'hbo'};
hb=job.run(hb);

% Get the link table with 
% bad short channels marked and filter for hbo, 
% copy over to hb probe
link_table = raw_data.probe.link(raw_data.probe.link.type == 850,:);
bad_idx = link_table.bad == 1;
hb.probe.link.bad(:) = 0;
hb.probe.link.bad(bad_idx) = 1;

% Calculate connectivity.
% ConnStats are sFCStats objects containing correlation values (all-to-all)
% This is the pearson correlation. It can also do coherence models.
job = nirs.modules.Connectivity;
job.AddShortSepRegressors = true;
ConnStats = job.run(hb);
Graph = ConnStats.graph('Z:hbo','p<0.005');

% Figures
% ConnStats.draw('R',[-1 1],'p<0.05');
% fname = 'ConnStats_' + subject_name + '.png';
% %saveas(gcf,fname)
% 
% figure('Position', get(0, 'Screensize'));
% Graph.pagerank.draw;
% fname = 'Graph_' + subject_name + '.png';
% %saveas(gcf,fname)
% 
% figure('Position', get(0, 'Screensize'));
% heatmap(ConnStats.R, 'Colormap', turbo);
% fname = 'R_heatmap_' + subject_name + '.png';
% %saveas(gcf,fname)

end