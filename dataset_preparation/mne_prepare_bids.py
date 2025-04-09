import os 
import matplotlib.pyplot as plt
import mne
import mne_nirs
import mne_bids


if __name__ == "__main__":

    # Should we plot the raw data for each subject?
    plot_subjects = False

    # Where are our files kept?
    original_folder = "../data/original_nirx_data"

    # Define a task name and a directory where to save the data to
    task = "rest"
    sessions = ["left", "right"]
    bids_root = "../data/park_move_rs_fnirs_bids"
    if not os.path.exists(bids_root):
        os.mkdir(bids_root)
    
    # Any subjects to skip?
    skip_subjects = []

    # Read the files and go through each subject
    nirx_folders = [f.path for f in os.scandir(original_folder) if f.is_dir()]
    for folder in nirx_folders:
        for session in sessions:

            # Get the ID of participant
            subject_id = os.path.basename(folder)
            print(f"Processing folder {folder} for subject {subject_id}")
            if subject_id in skip_subjects:
                print(f"Skipping {subject_id}")
                continue

            # Load the NIRS data for this subject
            folder_left = f"{folder}/{session}"
            raw_intensity = mne.io.read_raw_nirx(folder_left)

            # Plot
            if plot_subjects:
                no_channels = len(raw_intensity.ch_names)
                raw_intensity.plot(n_channels=no_channels, duration=500,
                                show_scrollbars=True, show_scalebars=True)
                raw_od = mne.preprocessing.nirs.optical_density(raw_intensity)
                raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)
                raw_haemo.plot(n_channels=no_channels, duration=500,
                            show_scrollbars=True, show_scalebars=True,
                            scalings=dict(hbo='1e-4', hbr='1e-4'))
                plt.show()

            # We'll have to convert the subject ID to BIDS-compliant standard
            # (no dashes etc)
            subject_id = subject_id.replace("_", "")

            # First convert to SNIRF and re-read the data from SNIRF
            mne_nirs.io.write_raw_snirf(raw_intensity, 'TEMP_snirf.snirf')
            snirf_intensity = mne.io.read_raw_snirf('TEMP_snirf.snirf', optode_frame='mri')

            # Write to BIDS format
            bids_path = mne_bids.BIDSPath(subject=subject_id, session=session, task=task, root=bids_root)
            mne_bids.write_raw_bids(snirf_intensity, bids_path, overwrite=True)
