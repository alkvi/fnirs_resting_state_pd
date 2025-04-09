% Ran on MATLAB R2024a
% BRAPH2.0 branch develop, latest commit Thu Nov 7 2024 
% commit sha 53095632a82c562f86e44ff415a613ce24f00f99

%% Settings

side = "left";

atlas_file_loc = 'atlasviewer_' + side + '.xlsx';
results_file_group = 'data/braph_' + side + '_group_results.csv';
results_file_subject = 'data/braph_' + side + '_subject_results.csv';
hc_group = 'HC_' + side;
pd_group = 'PD_' + side;

%% Load BrainAtlas

atlas_file = char(strcat('data', '/', 'braph_data', '/', atlas_file_loc));
atlas = readtable(atlas_file);
atlas_ch = atlas.Var1;

im_ba = ImporterBrainAtlasXLS( ...
    'FILE', atlas_file, ...
    'WAITBAR', true ...
    );

ba = im_ba.get('BA');

%% Load Groups

group_dir = char(strcat('data', filesep, 'braph_data', filesep, hc_group));
im_gr1 = ImporterGroupSubjectCON_XLS( ...
    'DIRECTORY', group_dir, ...
    'BA', ba, ...
    'WAITBAR', true ...
    );

gr1 = im_gr1.get('GR');

group_dir = char(strcat('data', filesep, 'braph_data', filesep, pd_group));
im_gr2 = ImporterGroupSubjectCON_XLS( ...
    'DIRECTORY', group_dir, ...
    'BA', ba, ...
    'WAITBAR', true ...
    );

gr2 = im_gr2.get('GR');

%% Group analysis

% Analyze Group 1 % Group 1 Analysis
a_WU1 = AnalyzeEnsemble_CON_WU('GR', gr1); 

% Analyze Group 2 % Group 2 Analysis
a_WU2 = AnalyzeEnsemble_CON_WU('GR', gr2, 'TEMPLATE', a_WU1); 

%% Run comparison

 % Compare Groups % Group Comparison
c_WU = CompareEnsemble('P', 1000, 'A1', a_WU1, 'A2', a_WU2, 'MEMORIZE', true);

ec_WU_diff = c_WU.get('COMPARISON', 'EigenVectorCentrality').get('DIFF');
ec_WU_p1 = c_WU.get('COMPARISON', 'EigenVectorCentrality').get('P1');
ec_WU_p2 = c_WU.get('COMPARISON', 'EigenVectorCentrality').get('P2');
ec_WU_cil = c_WU.get('COMPARISON', 'EigenVectorCentrality').get('CIL');
ec_WU_ciu = c_WU.get('COMPARISON', 'EigenVectorCentrality').get('CIU');

%% Save results

% Extract results
ec_WU_diff = ec_WU_diff{1};
ec_WU_p1 = ec_WU_p1{1};
ec_WU_p2 = ec_WU_p2{1};
ec_WU_cil = ec_WU_cil{1};
ec_WU_ciu = ec_WU_ciu{1};

% Also get the group averages behind the diffs
group_results_gr1 = a_WU1.get('ME_DICT').get('IT_LIST');
group_results_gr1 = group_results_gr1{1};
group_results_gr1 = group_results_gr1.get('M');
group_results_gr1 = group_results_gr1{1};
group_results_gr2 = a_WU2.get('ME_DICT').get('IT_LIST');
group_results_gr2 = group_results_gr2{1};
group_results_gr2 = group_results_gr2.get('M');
group_results_gr2 = group_results_gr2{1};

% Set up a table with results
result_table = table();
result_table.ch = atlas_ch;
result_table.hc_ec = group_results_gr1;
result_table.pd_ec = group_results_gr2;
result_table.diff = ec_WU_diff;
result_table.p1 = ec_WU_p1;
result_table.p2 = ec_WU_p2;
result_table.cil = ec_WU_cil;
result_table.ciu = ec_WU_ciu;

% Write results
writetable(result_table, results_file_group)

%% Subject results

% First get the subject level graphs from calculated group results
subject_results_gr1 = a_WU1.get('G_DICT').get('IT_LIST');
subject_results_gr2 = a_WU2.get('G_DICT').get('IT_LIST');
subject_list_gr1 = gr1.get('SUB_DICT').get('IT_LIST');
subject_list_gr2 = gr2.get('SUB_DICT').get('IT_LIST');

% Initialize an empty table and add channel list
subj_table = table();
subj_table.ch = atlas_ch;

% Go through each graph and extract the measurement dictionary
% and then the measure (we only have EigenvectorCentrality)
for subj_idx=1:length(subject_results_gr1)
    subject_measures = subject_results_gr1{subj_idx}.get('M_DICT').get('IT_LIST');
    subject_id = subject_list_gr1{subj_idx}.get('ID');
    subject_id = split(subject_id,'_');
    subject_id = subject_id{2};
    subject_ec = subject_measures{1}.get('M');
    subject_ec = subject_ec{1};

    % Add subject_ec as a column in the table with subject_id as the column name
    subj_table.(subject_id) = subject_ec;
end

for subj_idx=1:length(subject_results_gr2)
    subject_measures = subject_results_gr2{subj_idx}.get('M_DICT').get('IT_LIST');
    subject_id = subject_list_gr2{subj_idx}.get('ID');
    subject_id = split(subject_id,'_');
    subject_id = subject_id{2};
    subject_ec = subject_measures{1}.get('M');
    subject_ec = subject_ec{1};

    % Add subject_ec as a column in the table with subject_id as the column name
    subj_table.(subject_id) = subject_ec;
end

% Save the table
writetable(subj_table, results_file_subject);
