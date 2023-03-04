from external_imports import *
from image_process import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()
from scipy import stats

def average_recon_jv(path, voc_rad0):
    # Lists to hold data
    PLQEs = []
    fluxes = []
    num_suns = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad0/(sci.k*298))-1
    filenames = find_npy(f"{path}/PLQE_oc") 

    for filename in filenames:
        # Get the parameters
        num_sun = float(filename.split('_')[1])
        num_suns.append(num_sun)
        flux = float(filename.split('_')[2])
        fluxes.append(flux)

        PLQE =  np.mean(np.load(f"{path}/PLQE_oc/{filename}"))
        PLQEs.append(PLQE)

   
    bias = filename.split('_')[0] # Assuming constant bias at all points
    PLQEs = np.array(PLQEs)
    fluxes = np.array(fluxes)
    num_suns = np.array(num_suns)

    # calculate voc_rad intensity dependant
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    QFLSs =  voc_rads +  (sci.k*298/sci.e)*np.log(PLQEs)
    Js = (1-num_suns)

    #### sort the lists before returning, low to high
    QFLSs = np.array([i for _,i in sorted(zip(fluxes, QFLSs))])
    Js = np.array([i for _,i in sorted(zip(fluxes, Js))])
    num_suns.sort()
    fluxes.sort()


    return QFLSs, Js, bias, fluxes, num_suns


def intsweep_QFLS(path, voc_rad0, bias):
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad0/(sci.k*298))-1
    filenames = find_npy(f"{path}/PLQE_{bias}") 

    # save in the same place as the PLQE plots
    if not os.path.isdir(f"{path}/QFLS_int_{bias}"):
        os.makedirs(f"{path}/QFLS_int_{bias}")

    def process_one_file(filename):
        num_sun = float(filename.split('_')[1])
        flux = float(filename.split('_')[2])

        PLQE = np.load(f"{path}/PLQE_{bias}/{filename}")
        voc_rad = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_sun)+1)
        QFLS =  voc_rad +  (sci.k*298/sci.e)*np.log(PLQE)

        filename = f"{bias.upper()}_{num_sun}_{flux}_"
        np.save(f"{path}/QFLS_int_{bias}/{filename}", QFLS)
    Parallel(n_jobs=num_cores)(delayed(process_one_file)(filename) for filename in filenames)

def oc_sc_1sun(path, flux_1sun, device_area=0.3087):
    for i in ['oc', 'sc']:
        if not os.path.isdir(f"{path}/PLQE_oc_sc_1sun/{i}"):
            os.makedirs(f"{path}/PLQE_oc_sc_1sun/{i}")

    # For Voc_rad calculatuion:
    filenames_oc = find_npy(f"{path}/PLQE_oc") 
    fluxes = []
    for filename_oc in filenames_oc:
        flux = float(filename_oc.split('_')[2])
        fluxes.append(flux)
    fluxes.sort() # To make sure it's in some proper order

    sign_at_start = (fluxes[0] - flux_1sun)/abs(fluxes[0] - flux_1sun)
    for i, flux in enumerate(fluxes):
        if (flux - flux_1sun)/abs(flux - flux_1sun) != sign_at_start:
            flux_above = flux
            flux_bellow = fluxes[i-1]
            flux_above_index = i
            break

    for filename_oc in filenames_oc:
        filename_sc = 'SC_' + '_'.join(filename_oc.split('_')[1:])
        flux = float(filename_oc.split('_')[2])
        if flux == flux_above:
            PLQE_oc_above = np.load(f"{path}/PLQE_oc/{filename_oc}")
            PLQE_sc_above = np.load(f"{path}/PLQE_sc/{filename_sc}")
        elif flux == flux_bellow:
            PLQE_oc_bellow = np.load(f"{path}/PLQE_oc/{filename_oc}")
            PLQE_sc_bellow = np.load(f"{path}/PLQE_sc/{filename_sc}")
    
    # oc array interpolate:
    m_ocarr = (PLQE_oc_bellow - PLQE_oc_above)/ (flux_bellow-flux_above)
    c_ocarr = PLQE_oc_bellow - m_ocarr*flux_bellow
    oc_1sun = m_ocarr*flux_1sun + c_ocarr

    np.save(f"{path}/PLQE_oc_sc_1sun/oc/OC_1_{flux_1sun}", oc_1sun)
    # sc array interpolate:
    m_scarr = (PLQE_sc_bellow - PLQE_sc_above)/ (flux_bellow-flux_above)
    c_scarr = PLQE_sc_bellow - m_scarr*flux_bellow
    sc_1sun = m_scarr*flux_1sun + c_scarr
    np.save(f"{path}/PLQE_oc_sc_1sun/sc/SC_1_{flux_1sun}", sc_1sun)

    ####### get the sm stuff##########
    # OC
    vlist_temp = []
    jlist_temp = [] 
    with open(f"{path}/PLQE_oc/source_meter.csv" 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            vlist_temp.append(float(row[0]))
    Voc = inter([flux_above, flux_bellow], 
                [sorted(vlist_temp)[flux_above_index], sorted(vlist_temp)[flux_above_index-1]])(flux_1sun)
    with open(f"{path}/PLQE_sc/source_meter.csv" 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            jlist_temp.append(float(row[0]) * (1e3/device_area))
    Jsc = inter([flux_above, flux_bellow], 
                [sorted(jlist_temp)[flux_above_index], sorted(jlist_temp)[flux_above_index-1]])(flux_1sun)
    
    with open(f"{path}/PLQE_oc_sc_1sun/JV",'w') as file:
        writer = csv.writer()
        writer.writerows([['Voltage measured at oc (V)', Voc],
                          ['Current measured at sc (mAcm^-2)', Jsc]]) 

def intsweep_coleff(path):
    if not os.path.isdir(f"{path}/pseudo_col_eff"):
        os.makedirs(f"{path}/pseudo_col_eff")
    
    oc_filenames = find_npy(f"{path}/PLQE_oc")   
    for oc_filename in oc_filenames:
        sc_filename = 'SC_' + '_'.join(oc_filename.split('_')[1:])
        
        num_suns = float(oc_filename.split('_')[1])
        flux = float(oc_filename.split('_')[2])
        
        savename = f"NA_{num_suns}_{flux}_"

        #print(f"{oc_path}/{oc_file}")
        oc_image = np.load(f"{path}/PLQE_oc/{oc_filename}")
        sc_image = np.load(f"{path}/PLQE_sc/{sc_filename}")
        psudo_col_eff = (oc_image-sc_image)/oc_image
        
        np.save(f"{path}/pseudo_col_eff/{savename}",psudo_col_eff)

def get_idelity_map(path, sunmin, sunmax):

    lnmin = np.log(sunmin)
    lnmax = np.log(sunmax)

    filenames = find_npy(f"{path}/QFLS_int_oc")
    num_suns = np.array([float(i.split('_')[1]) for i in filenames])
    ln_curr =  np.log(num_suns)
    QFLS_arrs = [np.load(f"{path}/QFLS_int_oc/{filename}") for filename in filenames]

    n_id = np.zeros(QFLS_arrs[0].shape)
    # Must save so that we can memmap arr:
    np.save('temp', n_id)
    def process_row(i):
        # Opens arrat:
        arr = np.load('temp.npy', mmap_mode='r+')
        for j in range(QFLS_arrs[0].shape[1]):
            QFLS_pp = []
            for QFLS_arr in QFLS_arrs:
                QFLS_pp.append(QFLS_arr[i,j])
            
            QFLS_pp = [x for _,x in sorted(zip(ln_curr,QFLS_pp))]
            ln_curr.sort()

            minindex = np.argmin(abs(ln_curr-lnmin))
            maxindex = np.argmin(abs(ln_curr-lnmax))

            m,c = stats.linregress(ln_curr[minindex:maxindex], QFLS_pp[minindex:maxindex])[0:2]
            plt.scatter(ln_curr[minindex:maxindex], QFLS_pp[minindex:maxindex])
            plt.plot(ln_curr[minindex:maxindex], m*ln_curr[minindex:maxindex] + c)
               
            arr[i,j] = sci.e*m/(sci.k*293)
    Parallel(n_jobs=num_cores, verbose = 0)(delayed(process_row)(i) for i in range(QFLS_arrs[0].shape[0]))
    n_id_final = np.load('temp.npy')
    # Delete temp: 
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