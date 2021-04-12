from external_imports import *
from image_process import *

def vsweep_curr_map(path, savepath, cell_area=0.25):
    if not os.path.isdir(f"{savepath}\\vsweep_currmaps"):
        os.makedirs(f"{savepath}\\vsweep_currmaps")
    
    # Lists to hold data
    Vs = []
    Js = []
    with open(f"{path}\\vsweep\\source_meter.csv", 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            Vs.append(float(row[0]))
            Js.append(float(row[1])*1000/cell_area)

        
    filenames = find_npy(f"{path}\\vsweep") 

    # Loads sc inage:
    sc_filename = None
    for filename in filenames:
        if float(filename.split('_')[0]) == 0:
            sc_filename = filename
    if not sc_filename:
        raise FileNotFoundError("Couldn't find short circuit image!!")
    
    sc_im = np.load(f"{path}\\vsweep\\{sc_filename}")

    for filename in filenames:
        im =  np.mean(np.load(f"{path}\\vsweep\\{filename}"))
        im_volt = float(filename.split('_')[0])
        im_curr = Js[np.where(Vs == im_volt)[0][0]]

        curr_map = (im-sc_im) * im_curr / np.mean((im-sc_im))
        # save with voltages as the ID
        np.save(f"{savepath}\\vsweep_currmaps\\{im_volt}", curr_map)

    return Vs, Js
