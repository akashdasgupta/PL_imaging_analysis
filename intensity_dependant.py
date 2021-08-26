from external_imports import *
from image_process import *
from multiprocessing.pool import ThreadPool

def single_pix_recon_jv(path, row, col, voc_rad0):
    # Lists to hold data
    PLQEs = []
    fluxes = []
    num_suns = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad0/(sci.k*298))-1
    filenames = find_npy(f"{path}\\PLQE_oc") 

    for filename in filenames:
        # Get the parameters
        num_sun = float(filename.split('_')[1])
        num_suns.append(num_sun)
        flux = float(filename.split('_')[2])
        fluxes.append(flux)

        ave_rad = 3 # Good idea to average over a bunch of pixels
        PLQE =  np.mean(np.load(f"{path}\\PLQE_oc\\{filename}")[row-ave_rad:row+ave_rad,
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
    filenames = find_npy(f"{path}\\PLQE_oc") 

    for filename in filenames:
        # Get the parameters
        num_sun = float(filename.split('_')[1])
        num_suns.append(num_sun)
        flux = float(filename.split('_')[2])
        fluxes.append(flux)

        ave_rad = 3 # Good idea to average over a bunch of pixels
        PLQE =  np.mean(np.load(f"{path}\\PLQE_oc\\{filename}"))
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


def intsweep_QFLS(path, voc_rad0):
    PLQEs = []
    fluxes = []
    num_suns = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad0/(sci.k*298))-1
    filenames = find_npy(f"{path}\\PLQE_oc") 

    # save in the same place as the PLQE plots
    if not os.path.isdir(f"{path}\\QFLS_int_oc"):
        os.makedirs(f"{path}\\QFLS_int_oc")

    for filename in filenames:
        num_sun = float(filename.split('_')[1])
        num_suns.append(num_sun)
        flux = float(filename.split('_')[2])
        fluxes.append(flux)

        PLQE = np.load(f"{path}\\{filename}")
        PLQEs.append(np.mean(PLQE))
        voc_rad = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_sun)+1)
        QFLS =  voc_rad +  (sci.k*298/sci.e)*np.log(PLQE)

        filename = f"OC_{num_sun}_{flux}_"
        np.save(f"{path}\\QFLS_int_oc\\{filename}", QFLS)
    
    PLQEs = np.array(PLQEs)
    fluxes = np.array(fluxes)
    num_suns = np.array(num_suns)

    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)
    QFLSs = voc_rads +  (sci.k*298/sci.e)*np.log(PLQEs)
    Js = 1- num_suns
    return QFLSs, Js, flux, num_suns

def oc_sc_1sun(path, flux_1sun, device_area=0.3087):
    for i in ['oc', 'sc']:
        if not os.path.isdir(f"{path}\\PLQE_oc_sc_1sun\\{i}"):
            os.makedirs(f"{path}\\PLQE_oc_sc_1sun\\{i}")

    # For Voc_rad calculatuion:
    filenames_oc = find_npy(f"{path}\\PLQE_oc") 
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
        filename_sc = '_'.join(['SC']+filename_oc.split('_')[1:])
        flux = float(filename_oc.split('_')[2])
        if flux == flux_above:
            PLQE_oc_above = np.load(f"{path}\\PLQE_oc\\{filename_oc}")
            PLQE_sc_above = np.load(f"{path}\\PLQE_sc\\{filename_sc}")
        elif flux == flux_bellow:
            PLQE_oc_bellow = np.load(f"{path}\\PLQE_oc\\{filename_oc}")
            PLQE_sc_bellow = np.load(f"{path}\\PLQE_sc\\{filename_sc}")
    
    # oc array interpolate:
    m_ocarr = (PLQE_oc_bellow - PLQE_oc_above)/ (flux_bellow-flux_above)
    c_ocarr = PLQE_oc_bellow - m_ocarr*flux_bellow
    oc_1sun = m_ocarr*flux_1sun + c_ocarr
    np.save(f"{path}\\PLQE_oc_sc_1sun\\oc\\OC_1_{flux_1sun}", oc_1sun)
    # sc array interpolate:
    m_scarr = (PLQE_sc_bellow - PLQE_sc_above)/ (flux_bellow-flux_above)
    c_scarr = PLQE_sc_bellow - m_scarr*flux_bellow
    sc_1sun = m_scarr*flux_1sun + c_scarr
    np.save(f"{path}\\PLQE_oc_sc_1sun\\sc\\SC_1_{flux_1sun}", sc_1sun)

    ####### get the sm stuff##########
    # OC
    vlist_temp = []
    jlist_temp = [] 
    with open(f"{path}\\PLQE_oc\\source_meter.csv" 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            vlist_temp.append(float(row[0]))
    Voc = inter([flux_above, flux_bellow], 
                [sorted(vlist_temp)[flux_above_index], sorted(vlist_temp)[flux_above_index-1]])(flux_1sun)
    with open(f"{path}\\PLQE_sc\\source_meter.csv" 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            jlist_temp.append(float(row[0]) * (1e3/device_area))
    Jsc = inter([flux_above, flux_bellow], 
                [sorted(jlist_temp)[flux_above_index], sorted(jlist_temp)[flux_above_index-1]])(flux_1sun)
    
    with open(f"{path}\\PLQE_oc_sc_1sun\\JV",'w') as file:
        writer = csv.writer()
        writer.writerows([['Voltage measured at oc (V)', Voc],
                          ['Current measured at sc (mAcm^-2)', Jsc]]) 


    

        




# def jv_extrapolator(path, num_cores):
#     if not os.path.isdir(f"{path}\\V_inter"):
#         os.makedirs(f"{path}\\V_inter")
#     files = find_npy(f"{path}\\V")
    
#     arrmax, arrmin = (0, 10)
#     for file in files:
#         imarr = np.load(f"{path}\\V\\{file}")
#         if np.max(imarr) > arrmax:
#             arrmax = np.max(imarr)
#         if np.min(imarr) < arrmin:
#             arrmin = np.min(imarr)
#     vstep = (arrmax-arrmin)/100
#     vrange = np.arange(arrmin,arrmax+vstep,vstep)
    
#     arr_shape = imarr.shape
#     arrays = np.zeros((len(vrange), arr_shape[0], arr_shape[1]))
    
#     def single_pix(row,collength, path): 
#         for col in range(collength):
#             V = []
#             J = []
#             for file in files:
#                 J.append(float(file.split('_')[1]))
#                 vi = (np.load(f"{path}\\V\\{file}"))[row, col]
#                 V.append(vi)
#             f = inter(V,J,bounds_error=False)
#             arrays[:,row,col] = f(vrange)
    
#     Parallel(n_jobs=num_cores, backend='threading')(delayed(single_pix)(row, arr_shape[1], path) for row in range(arr_shape[0]))
    
#     for i,v in enumerate(vrange):
#         np.save(f"{path}\\V_inter\\{v}", arrays[i,:,:])