####################################
# UNVERIFIED that this works or produces real plots, DO NOT USE YET!!!
####################################

from external_imports import *
from image_process import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

def series_resistance(datapath, savepath, bandgap):
    # Load a low intensity OC image:

    
    vsweep_names = ['vsweep', "vsweep_f", "vsweep_b"]
    for vsweep_name in vsweep_names:
        if os.path.isdir(f"{datapath}\\{vsweep_name}"):
            break
    # Figure out MPP, Voc:
    Vs = []
    Js = []
    with open(f"{datapath}\\{vsweep_name}\\source_meter.csv",'r') as file:
        reader = csv.reader(file)
        for row in reader:
            Vs.append(float(row[0]))
            Js.append(float(row[1])*(1e3/0.3087))
    Js = np.array([i for _,i in sorted(zip(Vs, Js))])
    Vs = np.array(sorted(Vs))
   
    # np.argmin(Vs*Js) is near MPP, np.argmin(abs(Js)) gets close to Voc:
    int_mpp_V = Vs[np.argmin(Vs*Js)]
    Voc_approx = inter(Js, Vs)(0)
    midpoint_volt = (Voc_approx+int_mpp_V)/2
    
    diff_0 = 1e10
    for i in range(len(Vs)):
        diff_i = abs(Vs[i]-midpoint_volt)
        if diff_i < diff_0:
            diff_0 = diff_i
            V1 = Vs[i]
    
    print(V1)
    J1 = 1
    
    for filename in find_npy(f"{savepath}\\PLQE_{vsweep_name}"):
        V = float(filename.split('_')[0])
        if V == np.around(V1, 3):
            arr1 = np.load(f"{savepath}\\PLQE_{vsweep_name}\\{filename}")

    # OC loads:
    oc_filenames = find_npy(f"{savepath}\\PLQE_oc")
    num_suns = [float(i.split("_")[1]) for i in oc_filenames]
    oc_Js = [float(i.split("_")[2])*sci.e for i in oc_filenames]
    oc_filenames = [i for _,i in sorted(zip(num_suns, oc_filenames))]
    num_suns.sort()
    oc_Vs = []
    with open(f"{datapath}\\oc\\source_meter.csv",'r') as file:
        reader = csv.reader(file)
        for row in reader:
            oc_Vs.append(float(row[0]))
    oc_Vs.sort()
    focj = inter(oc_Vs, num_suns, bounds_error=False)
    
    oc_arrays = [np.load(f"{savepath}\\PLQE_oc\\{i}")*float(i.split('_')[1]) for i in oc_filenames]
    
    arr_shape = np.load(f"{savepath}\\PLQE_oc\\"+find_npy(f"{savepath}\\PLQE_oc")[0]).shape
    series_resistance_arr = np.zeros(arr_shape)
    np.save('temp', series_resistance_arr)
    
#     ## Figure out Jsc
#     sc_meas = []
#     with open(f"{datapath}\\sc\\source_meter.csv") as file:
#         reader = csv.reader(file)
#         for row in reader:
#             sc_meas.append(float(row[1])/(0.3087))
#     sc_meas.sort()
#     jsc = inter(num_suns, sc_meas)(1)
#     print(jsc)
    JG = j1sunf(1.6)*sci.e
    
    def process_row(i):
        # Opens array:
        arr = np.load('temp.npy', mmap_mode='r+')
        for j in range(arr_shape[1]):
            fV = inter([k[i,j] for k in oc_arrays], oc_Vs, bounds_error=False)
            V2 = fV(arr1[i,j])
            J2 = focj(V2)

            arr[i,j] = abs((V2-V1)/((J2-J1)*JG))
            
    
    Parallel(n_jobs=num_cores, verbose = 0)(delayed(process_row)(i) for i in range(arr_shape[0]))
    series_resistance_arr = np.load('temp.npy')
    # Delete temp: 
    os.remove('temp.npy')
    return series_resistance_arr