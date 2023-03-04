from external_imports import *
from PIL import Image

"""
These are functions that are not specific to the analysing data. Stuff like finding and opening files
"""


def path_process(path):
    """Finds directory structure of data:

    Returns dict with all substrates and repeats/pixels present, assuming data saved as:

    root path
    |_substrate_1
    | |_1
    | |_2
    | |_...
    |
    |_substrate_2
      |_...

    Args:
        path (str): root path to the data (with different substrates)

    Returns:
        db (dict): Dictonary of subdirs. Keys = substrates, elements = list of numbers in substrate
    """
    db = {}
    for root, dirs, _ in os.walk(path):
        for dir in dirs:
            if dir.lower() != "white":
                db[dir] = []
                for _, subdirs, _ in os.walk(f"{root}/{dir}"):
                    for subdir in subdirs:
                        try:
                            int(subdir)
                            db[dir].append(int(subdir))
                        except ValueError:
                            pass
                    break
        break
    return db


def find_npy(datapath):
    """Finds any file with the .npy extention in a path

    Args:
        datapath (str): root path to the data (with different substrates)

    Returns:
        raw_paths (list): list of filenames found ending in .npy
    """
    raw_paths = []
    for _, _, files in os.walk(datapath):
        for file in files:
            if file.endswith(".npy"):  # Only interested in np arrs
                raw_paths.append(file)
        break
    return raw_paths




