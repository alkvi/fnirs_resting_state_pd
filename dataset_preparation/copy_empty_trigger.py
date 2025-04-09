import os 
import shutil

if __name__ == "__main__":

    # Should we plot the raw data for each subject?
    plot_subjects = False

    # Where are our files kept?
    original_folder = "../data/original_nirx_data"

    # Define a task name and a directory where to save the data to
    task = "rest"
    bids_root = "../data/park_move_rs_fnirs_bids"
    if not os.path.exists(bids_root):
        os.mkdir(bids_root)

    # Any subjects to skip?
    skip_subjects = []

    # Read the files and go through each subject
    nirx_folders = [f.path for f in os.scandir(original_folder) if f.is_dir()]
    for folder in nirx_folders:

        # Get the ID of participant
        subject_id = os.path.basename(folder)
        print(f"Processing folder {folder} for subject {subject_id}")
        if subject_id in skip_subjects:
            print(f"Skipping {subject_id}")
            continue

        # Load the NIRS data for this subject
        folder_left = folder + "/left"
        orig_file = "../data/no-trig.tri"
        dest_file = folder_left + "/" + "no-trig.tri"
        print("Copying file to " + dest_file)
        shutil.copyfile(orig_file, dest_file)

        folder_left = folder + "/right"
        orig_file = "../data/no-trig.tri"
        dest_file = folder_left + "/" + "no-trig.tri"
        print("Copying file to " + dest_file)
        shutil.copyfile(orig_file, dest_file)