from external_imports import *
from PIL import Image

def path_process(path):
    db = {}  
    for root, dirs, _ in os.walk(path):
        for dir in dirs:
            if dir.lower() != 'white':
                db[dir] = []
                for _, subdirs, _ in os.walk(f"{root}\\{dir}"):
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
    raw_paths = []
    for _, _, files in os.walk(datapath):
        for file in files:
            if file.endswith(".tif"): # Only interested in tifs
                raw_paths.append(file)
        break
    return raw_paths

def find_npy(datapath):
    raw_paths = []
    for _, _, files in os.walk(datapath):
        for file in files:
            if file.endswith(".npy"): # Only interested in np arrs
                raw_paths.append(file)
        break
    return raw_paths


def get_cam_exposure(path):
    with open(path+'\\camera\\camera_setting_dump.txt') as file:
        for line in file:
            db = eval(line)
            db = eval(line)
            return float(db['ExposureTime'])

def image_intsweep_name_parser(path):
    tif_names = find_tif(path)
    repeat_db = {}

    for tif_name in tif_names:
        sm_output = tif_name.split('_')[0]
        led_v = tif_name.split('_')[1].split('=')[1]

        if f"{sm_output}_{led_v}" in repeat_db.keys():
            repeat_db[f"{sm_output}_{led_v}"].append(tif_name)
        else:
            repeat_db[f"{sm_output}_{led_v}"] = [tif_name]
    return repeat_db


def image_for_cropping(path):
    image_names = find_tif(path)

    nominal_vs = [i.split('_')[1].split('=')[1] for i in image_names]
    max_filename = image_names[np.argmax(nominal_vs)]

    im = np.array(Image.open(f"{path}\\{max_filename}"))
    im_scaled = im *(2**(8)/np.max(im)) # 8bit 
    image= im_scaled.astype(np.uint8)
    return image




