import os 
import shutil

if __name__ == "__main__":

    measure = "hbr"

    # Where all the data is stored
    snirf_source_folder = "data/subject_matrices"
    files = [ f.path for f in os.scandir(snirf_source_folder) if "_R.xlsx" in f.path]

    # Where we put data
    dest_folder_base = "data/braph_data"

    sides = ["left", "right"]
    groups = ["HC", "PD"]

    # Create side/measure/group folder
    for side in sides:
        for group in groups:
            
            dest_folder = dest_folder_base + "/" + group + "_" + side + "_" + measure
            print("Creating folder " + dest_folder)
            if not os.path.isdir(dest_folder):
                os.mkdir(dest_folder)

    # Copy over files to correct folder
    for file in files:
        print("Processing file " + file)
        file_base = os.path.basename(file)
        file_parts = file_base.split('_')
        subject = file_parts[1]
        side = file_parts[2]
        file_measure = file_parts[3]

        if file_measure != measure:
            continue

        group = "PD"
        if "FNP1" in subject:
            group = "HC"

        # Build the right path to send the file to
        dest_folder = dest_folder_base + "/" + group + "_" + side + "_" + measure
        dest_file = dest_folder + "/" + subject + "_" + side + "_" + group + "_" + measure + "_R.xlsx"

        # Copy file
        print("Copying file to " + dest_file)
        shutil.copyfile(file, dest_file)

    print("Done")