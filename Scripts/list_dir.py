from os import listdir
from os.path import isfile, join, isdir

def list_files(dir):
    onlyfiles = [ join(dir,f) for f in listdir(dir) if isfile(join(dir,f)) ]
    return onlyfiles

def list_subdirectories(dir):
    onlydirs = [ join(dir,f) for f in listdir(dir) if isdir(join(dir,f)) ]
    return onlydirs

if __name__ == "__main__":
    import sys
    for file in list_files(sys.argv[1]) :
        print(file)
