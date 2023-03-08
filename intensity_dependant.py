from external_imports import *
from image_process import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()
from scipy import stats

def intsweep_QFLS(path, voc_rad0, bias):
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad0/(sci.k*298))-1
    filenames = find_npy(f"{path}/PLQE_{bias}") 

    # save in the same place as the PLQE plots
    if not os.path.isdir(f"{path}/QFLS_{bias}"):
        os.makedirs(f"{path}/QFLS_{bias}")

    def process_one_file(filename):
        num_sun = float(filename.split('_')[1])
        flux = float(filename.split('_')[2])

        PLQE = np.load(f"{path}/PLQE_{bias}/{filename}")
        voc_rad = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_sun)+1)
        QFLS =  voc_rad +  (sci.k*298/sci.e)*np.log(PLQE)

        filename = f"{bias.upper()}_{num_sun}_{flux}_"
        np.save(f"{path}/QFLS_{bias}/{filename}", QFLS)
    Parallel(n_jobs=num_cores,backend='threading')(delayed(process_one_file)(filename) for filename in filenames)

    

def get_idelity_map(path, sunmin, sunmax, condition):

    lnmin = np.log(sunmin)
    lnmax = np.log(sunmax)

    filenames = find_npy(f"{path}/QFLS_{condition}")
    num_suns = np.array([float(i.split('_')[3]) for i in filenames])
    ln_curr =  np.log(num_suns)
    sorted_filenames = [x for _,x in sorted(zip(ln_curr,filenames))]
    ln_curr.sort()
    
    minindex = np.argmin(abs(ln_curr-lnmin))
    maxindex = np.argmin(abs(ln_curr-lnmax))
    ln_curr = np.array(ln_curr[minindex:maxindex])
    sorted_filenames = sorted_filenames[minindex:maxindex]
    
    QFLS_arrs = [np.load(f"{path}/QFLS_{condition}/{filename}") for filename in sorted_filenames]
    QFLS_arrs = np.array(QFLS_arrs)
    n_id = np.zeros(QFLS_arrs[0].shape)
    # Must save so that we can memmap arr:
    np.save('temp', n_id)
    
    indices = []
    for i in range(n_id.shape[0]):
        for j in range(n_id.shape[1]):
            indices.append((i,j))
    arr = np.load('temp.npy', mmap_mode='r+')  
    def process_row(index):
        i,j = index
        # Opens arrat:
        
        QFLS_pp = QFLS_arrs[:,i,j]
        m,c = stats.linregress(ln_curr, QFLS_pp)[0:2]
        arr[i,j] = sci.e*m/(sci.k*293)
            
    Parallel(n_jobs=num_cores, verbose = 0)(delayed(process_row)(index) for index in indices)
    del arr
    n_id_final = np.load('temp.npy')
#     # Delete temp: 
    os.remove('temp.npy')
    return n_id_final



def collect_meas_csv(rawpath, savepath, cell_area=0.3087):
    vocs = []
    jscs = []
    LEDs = []
    with open(f"{rawpath}/oc/source_meter.csv",'r') as file:
        reader = csv.reader(file)
        for row in reader:
            vocs.append(float(row[0]))
    with open(f"{rawpath}/oc/LED_power_supply.csv",'r') as file:
        reader = csv.reader(file)
        for row in reader:
            LEDs.append(float(row[0]))
    with open(f"{rawpath}/sc/source_meter.csv",'r') as file:
        reader = csv.reader(file)
        for row in reader:
            jscs.append(float(row[1])*(1e3/cell_area))
    ###assuming they all go low to high in the same way
    LEDs.sort()
    vocs.sort()
    jscs.sort()

    num_suns = [i.split('_')[1] for i in find_npy(f"{savepath}/PLQE_oc")]
    num_suns.sort()
    fluxes = [i.split('_')[2] for i in find_npy(f"{savepath}/PLQE_oc")]
    fluxes.sort()

    with open(f"{savepath}/Intensity_dependant_voc_Jsc.csv",'w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["LED voltage(V)", "Numsuns", "Flux (cm^-2)","Voc (V)", "Jsc (mAcm^-2)"])
        writer.writerows(zip(LEDs,num_suns,fluxes,vocs,jscs))