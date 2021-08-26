from external_imports import *
from image_process import *
from joblib import Parallel, delayed

def single_pix_recon_jv(path, row, col, voc_rad):
    # Lists to hold data
    PLQEs = []
    fluxes = []
    num_suns = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
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

    # calculate voc_rad intensity dependant
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    QFLSs =  voc_rads +  (sci.k*298/sci.e)*np.log(PLQEs)
    Js = (1-num_suns)
    
    return QFLSs, Js, bias, flux, num_suns


def average_recon_jv(path, voc_rad):
    # Lists to hold data
    PLQEs = []
    fluxes = []
    num_suns = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
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

    # calculate voc_rad intensity dependant
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    QFLSs =  voc_rads +  (sci.k*298/sci.e)*np.log(PLQEs)
    Js = (1-num_suns)
    
    return QFLSs, Js, bias, flux, num_suns


def array_jv(path, savepath, flux_1sun, voc_rad, white_mean_scaled):
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
        
    filenames_voc = find_npy(f"{path}\\oc") 

    if not os.path.isdir(f"{savepath}\\PLQE"):
        os.makedirs(f"{savepath}\\PLQE")
    if not os.path.isdir(f"{savepath}\\V"):
        os.makedirs(f"{savepath}\\V")

    for filename_voc in filenames_voc:
        # Asuming there is a sc file for every oc file, and it's saved by the same name format:
        flux = float(filename_voc.split('_')[1])
        exposure = float(filename_voc.split('_')[2])

        im_voc =  np.load(f"{path}\\oc\\{filename_voc}")/exposure
        # im_isc =  np.load(f"{path}\\sc\\{filename_isc}")/exposure
        
        white_ref = white_mean_scaled*flux 

        num_sun = flux/flux_1sun
        voc_rad = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_sun)+1)

        rr = (im_voc)/(white_ref)
        V =  voc_rad +  (sci.k*298/sci.e)*np.log(rr)
        # convention: [num suns]_[J]_[flux]_.npy
        filename = f"{num_sun}_{1-num_sun}_{flux}_"

        np.save(f"{savepath}\\PLQE\\PLQE_{filename}", rr)
        np.save(f"{savepath}\\V\\V{filename}", V)

def oc_sc_1sun(path, savepath, voc_rad, white_mean_scaled, led_1s):
    # For Voc_rad calculatuion:
    filenames_voc = find_npy(f"{path}\\oc") 
    flux_1sun = ledf(led_1s)

    for i, filename_voc in enumerate(filenames_voc):
        if i == len(filenames_voc):
            continue
        # Asuming there is a sc file for every oc file, and it's saved by the same name format:
        flux = float(filename_voc.split('_')[1])
        nextflux =  float(filenames_voc[i+1].split('_')[1])
        if (flux-flux_1sun)/abs(flux-flux_1sun) != (nextflux-flux_1sun)/abs(nextflux-flux_1sun):
            im_oc_before = np.load(f"{path}\\oc\\{filename_voc}")/float(filename_voc.split('_')[2])
            im_sc_before = np.load(f"{path}\\sc\\{'SC'+'_'.join(filename_voc.split('_')[1:])}")/float(filename_voc.split('_')[2])
           
            white_ref = white_mean_scaled*flux 
            PLQE_oc_before = (im_oc_before)/(white_ref)
            PLQE_sc_before = (im_sc_before)/(white_ref)
            
            im_after  = np.load(filenames_voc[i+1])/float(filenames_voc[i+1].split('_')[2])
            im_after  = np.load(filenames_voc[i+1])/float(filenames_voc[i+1].split('_')[2])
            break
            

    V =  voc_rad +  (sci.k*298/sci.e)*np.log(rr)
    # convention: [num suns]_[J]_[flux]_.npy
    filename = f"{num_sun}_{1-num_sun}_{flux}_"

    np.save(f"{savepath}\\PLQE\\PLQE_{filename}", rr)
    np.save(f"{savepath}\\V\\V{filename}", V)



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