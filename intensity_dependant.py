from external_imports import *
from image_process import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()
from scipy import stats

def single_pix_recon_jv(path, row, col, voc_rad0):
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

        ave_rad = 3 # Good idea to average over a bunch of pixels
        PLQE =  np.mean(np.load(f"{path}/PLQE_oc/{filename}")[row-ave_rad:row+ave_rad,
                                                                       col-ave_rad:col+ave_rad])
        PLQEs.append(PLQE)

   
    bias = filename.split('_')[0] # Assuming constant bias at all points
    PLQEs = np.array(PLQEs)
    fluxes = np.array(fluxes)
    num_suns = np.array(num_suns)

    # calculate voc_rad intensity dependant
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    QFLSs =  voc_rads +  (sci.k*298/sci.e)*np.log(PLQEs)
    Js = (1-num_suns)
    
    return QFLSs, Js, bias, flux, num_suns


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

def get_idelity_map(path):
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
            for kp in range(len(QFLS_arrs)-2):
                k=kp+1
                m,c = stats.linregress(ln_curr[0:k], QFLS_pp[0:k])[0:2]
                grad_vec = np.array([ln_curr[k+1]-ln_curr[k+1], QFLS_pp[k] +  m*(ln_curr[k+1]-ln_curr[k+1])])
                next_vector = np.array([ln_curr[k+1]-ln_curr[k+1], QFLS_pp[k+1]])
                dp = np.dot(grad_vec,next_vector)
                theta = np.arccos(dp/(np.sqrt(np.dot(grad_vec,grad_vec))*np.sqrt(np.dot(next_vector,next_vector))))
                if theta * 180/np.pi >= 2.5:
                    break
            arr[i,j] = sci.e*m/(sci.k*293)
    Parallel(n_jobs=num_cores, verbose = 0)(delayed(process_row)(i) for i in range(QFLS_arrs[0].shape[0]))
    n_id_final = np.load('temp.npy')
    # Delete temp: 
    os.remove('temp.npy')
    return n_id_final

def get_J0_map(path):
    filenames = find_npy(f"{path}/QFLS_int_oc")
    flux = np.array([float(i.split('_')[2]) for i in filenames])
    JGs =  flux * (sci.e*1e3) # in mA
    QFLS_arrs = [np.load(f"{path}/QFLS_int_oc/{filename}") for filename in filenames]
    
    plt.scatter(np.exp((sci.e*np.array([np.mean(i) for i in QFLS_arrs]))/(sci.k * 298)), JGs)
    plt.show()

    J0 = np.zeros(QFLS_arrs[0].shape)
    # Must save so that we can memmap arr:
    np.save('temp', J0)
    def process_row(i):
        # Opens arrat:
        arr = np.load('temp.npy', mmap_mode='r+')
        for j in range(QFLS_arrs[0].shape[1]):
            QFLS_pp = []
            for QFLS_arr in QFLS_arrs:
                QFLS_pp.append(QFLS_arr[i,j])
            
            QFLS_pp = np.array(QFLS_pp)
            exp_term = np.exp((sci.e*QFLS_pp)/(sci.k * 298))     
            exp_term = [i for _,i in sorted(zip(JGs,exp_term))]
            JGs.sort()
            
            for kp in range(len(JGs)-2):
                k=kp+1
                m,c = stats.linregress(exp_term[0:k], JGs[0:k])[0:2]
                grad_vec = np.array([exp_term[k+1]-exp_term[k+1], JGs[k] +  m*(exp_term[k+1]-exp_term[k+1])])
                next_vector = np.array([exp_term[k+1]-exp_term[k+1], JGs[k+1]])
                dp = np.dot(grad_vec,next_vector)
                theta = np.arccos(dp/(np.sqrt(np.dot(grad_vec,grad_vec))*np.sqrt(np.dot(next_vector,next_vector))))
                if theta * 180/np.pi >= 2.5:
                    break
            arr[i,j] = m


    Parallel(n_jobs=num_cores, verbose = 0)(delayed(process_row)(i) for i in range(QFLS_arrs[0].shape[0]))
    J0_final = np.load('temp.npy')
    # Delete temp: 
    os.remove('temp.npy')
    return J0_final

def get_Jsc_map(path):
    J0 = np.load(f"{path}/J0.npy")

    ### Get sc image at 1 sun
    filenames = find_npy(f"{path}/QFLS_int_sc")
    num_suns = np.array([float(i.split('_')[1]) for i in filenames])

    filenames = [i for _,i in sorted(zip(num_suns, filenames))]
    num_suns.sort()

    sign_0 = (num_suns[0] - 1)/abs(num_suns[0] - 1)
    for i,num_sun in enumerate(num_suns):
        if (num_sun-1)/abs(num_sun- 1) != sign_0:
            num_sun_above = num_sun
            arr_sc_above = np.load(f"{path}/QFLS_int_sc/{filenames[i]}")
            num_sun_bellow = num_suns[i-1]
            arr_sc_bellow = np.load(f"{path}/QFLS_int_sc/{filenames[i-1]}")
    
    m_arr = (arr_sc_above - arr_sc_bellow)/(num_sun_above-num_sun_bellow)
    c_arr = arr_sc_above - m_arr*num_sun_above
    QFLS_sc = m_arr + c_arr 
    ###

    JG = (float(filenames[0].split('_')[2]) /float(filenames[0].split('_')[1]))*(sci.e*1000)
    J_loss = J0 * np.exp((sci.e*QFLS_sc)/(sci.k * 298))
    J_sc = JG - J_loss
    return J_sc, J_loss