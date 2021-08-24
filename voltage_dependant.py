from external_imports import *
from image_process import *
from basic_funcs import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()


def vsweep_array_QFLS(path, savepath, voc_rad, white_mean_scaled, flux_1sun):
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
    filenames = find_npy(f"{path}\\vsweep") 

    if not os.path.isdir(f"{savepath}\\QFLS_vsweep"):
        os.makedirs(f"{savepath}\\QFLS_vsweep")

    for filename in filenames:
        # Asuming there is a sc file for every oc file, and it's saved by the same name format:
        flux = float(filename.split('_')[1])
        num_sun = flux / flux_1sun
        exposure = float(filename.split('_')[2])
        volt = filename.split('_')[0]

        im =  np.load(f"{path}\\vsweep\\{filename}")/exposure  
        im_volt = float(filename.split('_')[0].split('=')[1])
        white_ref = white_mean_scaled*flux 
        
        voc_rad = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_sun)+1)

        rr = (im)/(white_ref)
        V =  voc_rad +  (sci.k*298/sci.e)*np.log(rr)
        # convention: [num suns]_[J]_[flux]_.npy
        filename = f"{num_sun}_{flux}_"

        np.save(f"{savepath}\\QFLS_vsweep\\V{im_volt}_{filename}", V)

def vsweep_curr_map(path, savepath, cell_area=0.3087):
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

    if np.min(Js) <= 0 :
        sign = Js[0]/abs(Js[0])
        for i in range(len(Js)):
            if Js[i]/abs(Js[i]) != sign:
                oc_voltage_above = np.around(Vs[i], 3)
                J_above =  Js[i]
                oc_voltage_bellow = np.around(Vs[i-1], 3)
                J_bellow = Js[i-1]
                break
        for filename in filenames:
            if float(filename.split('_')[0].split('=')[1]) == oc_voltage_above:
                filename_above = filename
            elif float(filename.split('_')[0].split('=')[1]) == oc_voltage_bellow:
                filename_bellow = filename
            
        
        imarr2 = np.load(f"{path}\\vsweep\\{filename_bellow}") 
        imarr1 = np.load(f"{path}\\vsweep\\{filename_above}")
        oc_imarr = imarr1 - ((imarr2-imarr1)/(J_bellow-J_above))*J_above
    # elif Js[-1] <= 0:
    #     for i in range(len(Js)):
    #         if Js[len(Js)-1-i]>=0:
    #             oc_voltage_above = np.around(Vs[len(Js)-1-i], 3)
    #             oc_voltage_bellow = np.around(Vs[len(Js)-i], 3)
    #             break
    #     for filename in filenames:
    #         if float(filename.split('_')[0].split('=')[1]) == oc_voltage_above:
    #             filename_above = filename
    #         elif float(filename.split('_')[0].split('=')[1]) == oc_voltage_bellow:
    #             filename_bellow = filename
        
    #     oc_imarr = (np.load(f"{path}\\vsweep\\{filename_above}") + np.load(f"{path}\\vsweep\\{filename_bellow}"))/2
    
    else:
        oc_imarr = oc_image_maker(f"{path}\\oc",f"{path}\\vsweep")
    #oc_imarr = oc_image_maker(f"{path}\\oc",f"{path}\\vsweep")
    im_volts_sorted = []
    for filename in filenames:
        im_volts_sorted.append(float(filename.split('_')[0].split('=')[1]))
    im_volts_sorted.sort()

    for filename in filenames:
        im =  np.load(f"{path}\\vsweep\\{filename}")
        im_volt = float(filename.split('_')[0].split('=')[1])
        im_curr = Js[-im_volts_sorted.index(im_volt)]

        curr_map = (oc_imarr-im) / oc_imarr
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
        im_voc = m*(vsweep_flux-x[0])+y[0]


    elif  oc_flux_lists[-1] < vsweep_flux:
        x = []
        y = []
        for i in range(2):
            x.append(oc_flux_lists[len(oc_flux_lists)-i-1])
            for temp_filename in oc_file_list:
                if float(temp_filename.split('_')[1]) == oc_flux_lists[len(oc_flux_lists)-i-1]:
                    y.append(np.load(f"{oc_path}\\{temp_filename}")/ float(temp_filename.split('_')[2]))
        
        m = (y[0]-y[1])/ x[0]-x[1]
        im_voc = m*(vsweep_flux-x[0])+y[0]

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

def int_col_eff(path, savepath, flux_1sun):
    
    if not os.path.isdir(f"{savepath}\\pseudo_col_eff"):
        os.makedirs(f"{savepath}\\pseudo_col_eff")
    
    oc_path = f"{path}\\oc"
    sc_path = f"{path}\\sc"
    
    for _,_,files in os.walk(oc_path):
        for oc_file in files:
            #print(files)
            if oc_file.endswith('.npy'):
                sc_file = 'SC_' + '_'.join(oc_file.split('_')[1:])
                
                flux = float(oc_file.split('_')[1])
                num_sun = flux/flux_1sun
                filename = f"coleff_{num_sun}_{1-num_sun}_{flux}_"

                #print(f"{oc_path}\\{oc_file}")
                oc_image = np.load(f"{oc_path}\\{oc_file}")
                sc_image = np.load(f"{sc_path}\\{sc_file}")
                psudo_col_eff = (oc_image-sc_image)/oc_image
                
                
                np.save(f"{savepath}\\pseudo_col_eff\\{filename}",psudo_col_eff)


