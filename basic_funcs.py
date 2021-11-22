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


def find_tif(datapath):
    """Finds any file with the .tif extention in a path

    Args:
        datapath (str): root path to the data (with different substrates)

    Returns:
        raw_paths (list): list of filenames found ending in .tif
    """
    raw_paths = []
    for _, _, files in os.walk(datapath):
        for file in files:
            if file.endswith(".tif"):  # Only interested in tifs
                raw_paths.append(file)
        break
    return raw_paths


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


def get_cam_exposure(path):
    """Pulls camera exposure time from setting dump file

    Args:
        path (str): path where camera data is located (NOT where dump file is)

    Returns:
        (float) : exposure time
    """
    with open(path + "/camera/camera_setting_dump.txt") as file:
        for line in file:
            db = eval(line)
            db = eval(line)
            return float(db["ExposureTime"])


def image_intsweep_name_parser(path):
    """Returns dict of repeat .tif data

    Args:
        path (str): Where your data lives

    Returns:
        repeat_db (dict): key = bias conditions(electrical, optical), element = list of filenames that are repeats
    """
    tif_names = find_tif(path)
    repeat_db = {}

    for tif_name in tif_names:
        sm_output = tif_name.split("_")[0]
        led_v = tif_name.split("_")[1].split("=")[1]

        if f"{sm_output}_{led_v}" in repeat_db.keys():
            repeat_db[f"{sm_output}_{led_v}"].append(tif_name)
        else:
            repeat_db[f"{sm_output}_{led_v}"] = [tif_name]
    return repeat_db


def image_for_cropping(path):
    """Pulls up an image from the data so the interactive cropper has something to show

    Args:
        path (str): where the data lives

    Returns:
        image (numpy.ndarray): image (scaled to 8 bit) for cropping
    """
    image_names = find_tif(path)

    nominal_vs = [i.split("_")[1].split("=")[1] for i in image_names]
    max_filename = image_names[np.argmax(nominal_vs)]

    im = np.array(Image.open(f"{path}/{max_filename}"))
    im_scaled = im * (2 ** (8) / np.max(im))  # 8bit
    image = im_scaled.astype(np.uint8)
    return image
