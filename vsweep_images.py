from external_imports import *
from image_process import *
from basic_funcs import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()


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

    if Js[0] <= 0 :
        for i in range(len(Js)):
            if Js[i]>=0:
                oc_voltage_above = np.around(Vs[i], 3)
                oc_voltage_bellow = np.around(Vs[i-1], 3)
                break
        for filename in filenames:
            if float(filename.split('_')[0].split('=')[1]) == oc_voltage_above:
                filename_above = filename
            elif float(filename.split('_')[0].split('=')[1]) == oc_voltage_bellow:
                filename_bellow = filename
        
        oc_imarr = (np.load(f"{path}\\vsweep\\{filename_above}") + np.load(f"{path}\\vsweep\\{filename_bellow}"))/2
    
    elif Js[-1] <= 0:
        for i in range(len(Js)):
            if Js[len(Js)-1-i]>=0:
                oc_voltage_above = np.around(Vs[len(Js)-1-i], 3)
                oc_voltage_bellow = np.around(Vs[len(Js)-i], 3)
                break
        for filename in filenames:
            if float(filename.split('_')[0].split('=')[1]) == oc_voltage_above:
                filename_above = filename
            elif float(filename.split('_')[0].split('=')[1]) == oc_voltage_bellow:
                filename_bellow = filename
        
        oc_imarr = (np.load(f"{path}\\vsweep\\{filename_above}") + np.load(f"{path}\\vsweep\\{filename_bellow}"))/2
    
    else:
        oc_imarr = oc_image_maker(f"{path}\\oc",f"{path}\\vsweep")

    im_volts_sorted = []
    for filename in filenames:
        im_volts_sorted.append(float(filename.split('_')[0].split('=')[1]))
    im_volts_sorted.sort()

    for filename in filenames:
        im =  np.load(f"{path}\\vsweep\\{filename}")
        im_volt = float(filename.split('_')[0].split('=')[1])
        im_curr = Js[-im_volts_sorted.index(im_volt)]

        curr_map = (oc_imarr-im) * im_curr / np.mean((oc_imarr-im))
        # save with voltages as the ID
        np.save(f"{savepath}\\vsweep_currmaps\\{im_volt}", curr_map)


    return Vs, Js


def oc_image_maker(oc_path, vsweep_path):
    vsweep_exposure = float(find_npy(vsweep_path)[0].split('_')[2])
    vsweep_flux = float(find_npy(vsweep_path)[0].split('_')[1])

    oc_flux_lists = []
    oc_file_list = find_npy(oc_path)
    for filename in oc_file_list:
        oc_flux_lists.append(float(filename.split('_')[1]))
        oc_flux_lists.sort()
    #print(oc_flux_lists)
    
    #im_voc = None
    if oc_flux_lists[0] > vsweep_flux:
        x = []
        y = []
        for i in range(2):
            x.append(oc_flux_lists[i])
            for temp_filename in oc_file_list:
                if float(temp_filename.split('_')[1]) == oc_flux_lists[i]:
                    y.append(np.load(f"{oc_path}\\{temp_filename}")/ float(temp_filename.split('_')[2]))
        m = (y[0]-y[1])/ x[0]-x[1]
        im_voc = m*(-1*x[0])+y[0]


    elif  oc_flux_lists[-1] < vsweep_flux:
        x = []
        y = []
        for i in range(2):
            x.append(oc_flux_lists[len(oc_flux_lists)-i-1])
            for temp_filename in oc_file_list:
                if float(temp_filename.split('_')[1]) == oc_flux_lists[len(oc_flux_lists)-i-1]:
                    y.append(np.load(f"{oc_path}\\{temp_filename}")/ float(temp_filename.split('_')[2]))
        
        m = (y[0]-y[1])/ x[0]-x[1]
        im_voc = m*(-1*x[0])+y[0]

    else:
        for i, oc_flux in enumerate(oc_flux_lists):
            if oc_flux>vsweep_flux:
                break
        for filename in find_npy(oc_path):
            if float(filename.split('_')[1]) == oc_flux_lists[i]:
                im_above = np.load(f"{oc_path}\\{filename}") / float(filename.split('_')[2])
            elif float(filename.split('_')[1]) == oc_flux_lists[i-1]:
                im_bellow = np.load(f"{oc_path}\\{filename}") / float(filename.split('_')[2])
        im_voc = (im_above + im_bellow) / 2
    
    im_voc *= vsweep_exposure
    return im_voc


