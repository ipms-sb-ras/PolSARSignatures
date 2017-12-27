from list_dir import list_files, list_subdirectories
from plot_results import plot_3d_with_contour, plot_3d_and_contour
import os


def save_all_in_folder(foldername, column):
    files = list_files(foldername)
    print('Number of files {}'.format(len(files)))
    if len(files) == 0:
        print('Folder "{}" is empty.'.format(foldername))
        return None

    ext = '.avrg'
    for f in files:
        if f.endswith(ext):
            print(f)
            # plot_3d_with_contour(f, column)
            plot_3d_and_contour(f, column)


def save_all_in_subfolders(parentfolder, column):
    for root, dirs, files in os.walk(parentfolder):
        for dir in dirs:
            path = os.path.join(root,dir)
            print('Processing folder {}'.format(path))
            save_all_in_folder(path, column)

if __name__ == "__main__":
    import sys
    save_all_in_subfolders(sys.argv[1], int(sys.argv[2]))
