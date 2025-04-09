import os
import mne
import numpy as np
import pandas as pd
from itertools import compress
import matplotlib.pyplot as plt

import mne_nirs
from mne.preprocessing.nirs import optical_density
from mne_bids import BIDSPath, read_raw_bids, print_dir_tree, make_report


if __name__ == "__main__":

    # Where will results go
    result_file_name = "../data/bad_shorts.csv"

    # Specify BIDS root folder
    bids_root = "../data/park_move_rs_fnirs_bids"
    nirx_folder = "../data/original_nirx_data"

    plot_subject = True

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
    sessions = ["left", "right"]

    # Have we already written some results?
    data_exists = False
    existing_data = None
    if os.path.isfile(result_file_name):
        data_exists = True
        existing_data = pd.read_csv(result_file_name)
        print(existing_data)

    # Go through data for each subject
    subject_frames = []
    for subject in all_subjects:
        for session in sessions:

            # Does this session already exist?
            processed = (data_exists) and \
            (not existing_data[(existing_data['subject'] == subject) & (existing_data['session'] == session)].empty)
            if processed:
                print(f"Skipping processed session {subject} / {session}")
                continue

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
            ch_names = np.array(raw_intensity.ch_names)

            # Calculate overall SCI over entire length of signal
            od = optical_density(raw_intensity)
            overall_sci = mne.preprocessing.nirs.scalp_coupling_index(od)
            overall_sci = np.array(overall_sci)
            overall_sci_short = overall_sci[short_idx]

            # Get bad channels for subject/session
            bad_channels = list(compress(raw_intensity.ch_names, overall_sci < 0.7))
            od.info["bads"] = list(compress(od.ch_names, overall_sci < 0.7))
            bad_channels = np.unique([channel.replace(" 760", "").replace(" 850", "") for channel in np.unique(bad_channels)]).tolist()
            print(bad_channels)

            # MBLL
            raw_haemo = mne.preprocessing.nirs.beer_lambert_law(od, ppf=0.1)

            # Get only short channels
            raw_haemo_short = raw_haemo.copy()
            picks = mne.pick_types(raw_haemo_short.info, meg=False, fnirs=True)
            dists = mne.preprocessing.nirs.source_detector_distances(raw_haemo_short.info, picks=picks)
            raw_haemo_short = raw_haemo_short.pick(picks[dists < 0.01])

            # Load aux data
            # First find matching NIRx folder
            nirx_subj_folder = [f.path for f in os.scandir(nirx_folder) if (f.is_dir() and subject.replace('rs', '') in f.path)]
            nirx_subj_folder = f"{nirx_subj_folder[0]}/{session}"
            nirx_snirf_file = [f.path for f in os.scandir(nirx_subj_folder) if "snirf" in f.path][0]
            raw_nirx_snirf = mne.io.read_raw_snirf(nirx_snirf_file)
            aux_df = mne_nirs.io.read_snirf_aux_data(nirx_snirf_file, raw_nirx_snirf)

            # Short channel plots
            no_channels = len(raw_haemo_short.ch_names)
            raw_haemo_short.plot(n_channels=no_channels, duration=500, picks="hbo",
                        show_scrollbars=True, show_scalebars=True,
                        scalings=dict(hbo='1e-4', hbr='1e-4'))
            
            fig = raw_haemo_short.compute_psd().plot(average=False, amplitude=False, picks="data")
            fig.suptitle(f"PSD shorts", weight="bold", size="x-large")

            # Ask for bad channels to mark
            bad_short_string = input("Enter bad channels S#-D#, separate with comma (,):\n")
            plt.show()

            short_data = {
                'subject': [subject],
                'session': [session],
                'bad_shorts': [bad_short_string],
            }
            short_frame = pd.DataFrame(data=short_data)
            subject_frames.append(short_frame)

    # Concatenate all frames
    full_frame = pd.concat(subject_frames)
    print(full_frame)

    # Write to csv
    if data_exists:
        full_frame.to_csv(result_file_name, index=False, mode='a', header=False)
    else:
        full_frame.to_csv(result_file_name, index=False)

    