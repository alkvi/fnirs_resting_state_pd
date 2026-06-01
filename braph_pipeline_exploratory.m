% Ran on MATLAB R2024a
% BRAPH2.0 branch develop, latest commit Thu Nov 7 2024
% commit sha 53095632a82c562f86e44ff415a613ce24f00f99

%% Settings

rng(1,"twister");
side = "left";

atlas_file_loc = 'data_public/atlasviewer_' + side + '.xlsx';
results_file_global = 'data/braph_' + side + '_exploratory_results_global.csv';
results_file_nodal = 'data/braph_' + side + '_exploratory_results_nodal.csv';
results_file_subject = 'data/braph_' + side + '_exploratory_subject_results.csv';
hc_group = 'HC_' + side;
pd_group = 'PD_' + side;

%% Load BrainAtlas

atlas_file = char(atlas_file_loc);
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

% Global results in a 1x1 cell with 1x1 double
% Nodal results in a 1x1 cell with 42x1 double
% Binodal results in 1x1 cell with 42x42 double
global_measures = {
    'Assortativity', 'ClusteringAv', 'DegreeAv', 'Diameter', 'EccentricityAv', ...
    'GlobalEfficiencyAv', 'LocalEfficiencyAv', 'Modularity', 'PathLengthAv', ...
    'Radius', 'RichClub', 'SmallWorldness', 'StrengthAv', 'Transitivity'
};

nodal_measures = {
    'BetweennessCentrality', 'Clustering', 'CommunityStructure', 'CorePeriphery', ...
    'Degree', 'Eccentricity', 'EigenVectorCentrality', 'GlobalEfficiency', ...
    'KCorenessCentrality', 'LocalEfficiency', 'Participation', 'PathLength', ...
    'RCDeg', 'RCS', 'Richness', 'Strength', 'Triangles'
};

global_results = table();
nodal_results = table();

c_WU = CompareEnsemble('P', 1000, 'A1', a_WU1, 'A2', a_WU2, 'MEMORIZE', true);

for i = 1:length(global_measures)
    measure = global_measures{i};
    diff = c_WU.get('COMPARISON', measure).get('DIFF');
    p1 = c_WU.get('COMPARISON', measure).get('P1');
    p2 = c_WU.get('COMPARISON', measure).get('P2');
    cil = c_WU.get('COMPARISON', measure).get('CIL');
    ciu = c_WU.get('COMPARISON', measure).get('CIU');

    global_results = [global_results; table({measure}, diff{1}, p1{1}, p2{1}, cil{1}, ciu{1}, ...
        'VariableNames', {'Measure', 'Diff', 'P1', 'P2', 'CIL', 'CIU'})];
end
writetable(global_results, results_file_global);


for i = 1:length(nodal_measures)
    measure = nodal_measures{i};
    diff = c_WU.get('COMPARISON', measure).get('DIFF');
    diff = diff{1};
    p1 = c_WU.get('COMPARISON', measure).get('P1');
    p1 = p1{1};
    p2 = c_WU.get('COMPARISON', measure).get('P2');
    p2 = p2{1};
    cil = c_WU.get('COMPARISON', measure).get('CIL');
    cil = cil{1};
    ciu = c_WU.get('COMPARISON', measure).get('CIU');
    ciu = ciu{1};

    % Repeat the measure name and channel labels
    measure_col = repmat({measure}, length(diff), 1);
    ch_col = atlas_ch(:);

    % Append to nodal_results
    nodal_results = [nodal_results; table(measure_col, ch_col, diff, p1, p2, cil, ciu, ...
        'VariableNames', {'Measure', 'Channel', 'Diff', 'P1', 'P2', 'CIL', 'CIU'})];
end
writetable(nodal_results, results_file_nodal);

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
% and then the measure
for subj_idx=1:length(subject_results_gr1)
    subject_measures = subject_results_gr1{subj_idx}.get('M_DICT').get('IT_LIST');
    subject_id = subject_list_gr1{subj_idx}.get('ID');
    subject_id = split(subject_id,'_');
    subject_id = subject_id{2};

    % Loop through all measures (1x31 cell)
    for m_idx = 1:length(subject_measures)
        measure_obj = subject_measures{m_idx};

        % Get measure name (e.g. 'Assortativity', 'EigenvectorCentrality')
        measure_name = measure_obj.get('ID');  

        if ~ismember(measure_name, nodal_measures)
            continue;
        end

        % Get numeric value
        measure_value = measure_obj.get('M');
        measure_value = measure_value{1};

        % Create column name: subjectID_measureName
        col_name = [subject_id '_' measure_name];

        % Store in table
        subj_table.(col_name) = measure_value;
    end

end

for subj_idx=1:length(subject_results_gr2)
    subject_measures = subject_results_gr2{subj_idx}.get('M_DICT').get('IT_LIST');
    subject_id = subject_list_gr2{subj_idx}.get('ID');
    subject_id = split(subject_id,'_');
    subject_id = subject_id{2};

    % Loop through all measures (1x31 cell)
    for m_idx = 1:length(subject_measures)
        measure_obj = subject_measures{m_idx};

        % Get measure name (e.g. 'Assortativity', 'EigenvectorCentrality')
        measure_name = measure_obj.get('ID');

        if ~ismember(measure_name, nodal_measures)
            continue;
        end

        % Get numeric value
        measure_value = measure_obj.get('M');
        measure_value = measure_value{1};

        % Create column name: subjectID_measureName
        col_name = [subject_id '_' measure_name];

        % Store in table
        subj_table.(col_name) = measure_value;
    end

end

% Save the table
writetable(subj_table, results_file_subject);
