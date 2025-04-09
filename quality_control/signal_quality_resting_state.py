import os
import mne
import pickle
import json
import numpy as np
from itertools import compress
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()

import mne_nirs
from mne.preprocessing.nirs import optical_density
from mne_bids import BIDSPath, read_raw_bids, print_dir_tree, make_report
from collections import defaultdict


def plot_bar(raw_intensity, sci_data, short_idx, fig_folder, fig_name):
    channel_names = raw_intensity.ch_names
    channel_names = np.delete(channel_names, short_idx)
    sns.set(rc = {'figure.figsize':(15,10)})
    sns.set_color_codes("pastel")
    ax = sns.barplot(x=sci_data, y=channel_names, data=None, color="b")
    ax.axvline(x=0.7, color='r', linestyle='dashed')
    ax.set(xlim=(0, 1), ylabel="",
        xlabel="SCI value")
    sns.despine(left=True, bottom=True)
    plt.savefig(fig_folder + "/" + fig_name)
    plt.close()

if __name__ == "__main__":

    # Specify BIDS root folder
    bids_root = "../data/park_move_rs_fnirs_bids"
    #print(make_report(bids_root))

    plot_subject = False

    # We have 4 file types: events, channels, optodes, and nirs.
    datatype = 'nirs'
    bids_path = BIDSPath(root=bids_root, datatype=datatype)
    nirs_files = bids_path.match()

    # For visualization
    subjects_dir = str(mne.datasets.sample.data_path()) + '/subjects'
    mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir)
    mne.set_config('MNE_NIRS_FOLD_PATH', 'C:/Users/alkvis/OneDrive - Karolinska Institutet/Dokument/Software/fOLD-public/Supplementary')

    # Get all subjects in a sorted list
    all_subjects = [file.subject for file in nirs_files]
    all_subjects = sorted(list(set(all_subjects)))
    print("Subjects")
    print(all_subjects)

    # Prepare params
    task = 'rest'
    suffix = 'nirs'
    session = 'left'

    # Save data in dicts for each protocol/condition
    all_sci_per_participant_long = {}
    all_sci_per_participant_short = {}
    all_sci_bad_channels_per_participant = {}
    sci_array = []
    bad_ch_array = []
    power_array = []

    # Go through data for each subject
    for subject in all_subjects:

        #if subject not in ["FNP1011rs", "FNP1012rs"]:
        #    continue

        # Load
        bids_path = BIDSPath(subject=subject, task=task, session=session,
                            suffix=suffix, datatype=datatype, root=bids_root)
        print("Using BIDS file path..")
        print(bids_path)
        if not os.path.exists(bids_path):
            print("No SNIRF file, skip: " + str(bids_path))
            continue
        raw_intensity = read_raw_bids(bids_path=bids_path, verbose=False)

        # Which channels are short channels? Have 16 (8 for 850, 8 for 760)
        short_idx = mne.preprocessing.nirs.short_channels(raw_intensity.info)
        short_idx = np.where(short_idx)[0]
        long_ch_names = np.array(raw_intensity.ch_names)
        long_ch_names = np.delete(long_ch_names, short_idx)

        # Calculate overall SCI over entire length of signal
        od = optical_density(raw_intensity)
        overall_sci = mne.preprocessing.nirs.scalp_coupling_index(od)
        overall_sci = np.array(overall_sci)
        overall_sci_short = overall_sci[short_idx]
        overall_sci_long = np.delete(overall_sci, short_idx, axis=0)
        sci_array.append(overall_sci)

        # Calculate peak power
        _, scores, _ = mne_nirs.preprocessing.peak_power(od, time_window=10)
        peak_power = np.mean(scores, axis=1)
        power_array.append(peak_power)
        
        # Add to dict
        subj_key = f"{subject}"
        all_sci_per_participant_long[subj_key] = [overall_sci_long]
        all_sci_per_participant_short[subj_key] = [overall_sci_short]

        # Get bad channels for subject/session
        bad_channels = list(compress(raw_intensity.ch_names, overall_sci < 0.7))
        bad_channels = np.unique([channel.replace(" 760", "").replace(" 850", "") for channel in np.unique(bad_channels)]).tolist()
        if subj_key not in all_sci_bad_channels_per_participant:
            all_sci_bad_channels_per_participant[subj_key] = []
        if len(bad_channels) > 0:
            all_sci_bad_channels_per_participant[subj_key] = bad_channels

        # Store logical list of bad channel idx in array
        bad_ch_array.append(overall_sci < 0.7)   

        if plot_subject:

            _, scores, times = mne_nirs.preprocessing.scalp_coupling_index_windowed(od, time_window=60)
            mne_nirs.visualisation.plot_timechannel_quality_metric(
                od,
                scores,
                times,
                threshold=0.7,
                title="Scalp Coupling Index " "Quality Evaluation",
            )
            od, scores, times = mne_nirs.preprocessing.peak_power(od, time_window=10)
            mne_nirs.visualisation.plot_timechannel_quality_metric(
                od, scores, times, threshold=0.1, title="Peak Power Quality Evaluation"
            )
            fig = mne_nirs.visualisation.plot_nirs_source_detector(
                overall_sci,
                raw_intensity.info,
                surfaces="brain",
                subject="fsaverage",
                subjects_dir=subjects_dir,
                trans="fsaverage",
            )
            plt.show()
            input("PRESS ENTER TO CONTINUE")

        # Plot overall SCI
        fig_folder = "figures"
        fig_name = f"SCI_overall_{bids_path.basename}_{session}.png"
        plot_bar(raw_intensity, overall_sci_long, short_idx, fig_folder, fig_name)

sci_array = np.array(sci_array)
sci_array_avg = np.mean(sci_array, axis=0)
bad_ch_array = np.array(bad_ch_array)
percent_bad_ch = np.sum(bad_ch_array, axis=0) / len(all_subjects) * 100

power_array = np.array(power_array)
power_array_avg = np.mean(power_array, axis=0)

view_map = {
    'left-lat': np.arange(1, 56),
    'left-frontal': np.arange(1, 56),
    'right-frontal': np.arange(1, 56),
    'right-lat': np.arange(1, 56),
}
fig_montage = mne_nirs.visualisation.plot_3d_montage(
    raw_intensity.info,
    view_map=view_map,
    subjects_dir=subjects_dir,
    src_det_names=None,
    ch_names=defaultdict(lambda: '')
)

# SCI, average
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    sci_array_avg,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
# Bad channels, number of
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    percent_bad_ch,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
# Peak power, average
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    power_array_avg,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
plt.show()
input("PRESS ENTER TO CONTINUE")

print("SCI per participant")
for key in all_sci_per_participant_long:
    print(f"\nKey: {key}")
    sci_for_subj_long = all_sci_per_participant_long[key][0]
    sci_for_subj_short = all_sci_per_participant_short[key][0]
    print("\nMean SCI value")
    print(np.mean(sci_for_subj_long))
    print(np.mean(sci_for_subj_short))

    print("\nMean SCI values, SD")
    print(np.std(sci_for_subj_long))
    print(np.std(sci_for_subj_short))

    sci_vals_bad_long = sci_for_subj_long[sci_for_subj_long<=0.7]
    sci_vals_bad_short = sci_for_subj_short[sci_for_subj_short<=0.7]
    print("\nCount of bad SCI values (<0.7)")
    print(len(sci_vals_bad_long))
    print(len(sci_vals_bad_short))