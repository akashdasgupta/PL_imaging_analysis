from external_imports import *
from image_process import *

def single_pixel_jv(path, row, col, flux_1sun, voc_rad, white_mean_scaled):
    # Lists to hold data
    dp_oc = []
    dp_sc = []
    white_refs = []
    fluxes = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
        
    filenames_voc = find_npy(f"{path}\\oc")    


    for filename_voc in filenames_voc:
        # Asuming there is a sc file for every oc file, and it's saved by the same name format:
        filename_isc = 'SC_'+'_'.join(filename_voc.split('_')[1:])
        flux = float(filename_voc.split('_')[1])
        exposure = float(filename_voc.split('_')[2])
        fluxes.append(flux)

        ave_rad = 3 # Good idea to average over a bunch of pixels

        datapoint_voc =  np.mean(np.load(f"{path}\\oc\\{filename_voc}")[row-ave_rad:row+ave_rad,
                                                                       col-ave_rad:col+ave_rad])
        datapoint_isc =  np.mean(np.load(f"{path}\\sc\\{filename_isc}")[row-ave_rad:row+ave_rad,
                                                                       col-ave_rad:col+ave_rad])
        
        dp_oc.append(datapoint_voc/exposure)
        dp_sc.append(datapoint_isc/exposure)
        white_refs.append(white_mean_scaled*flux)

    dp_oc = np.array(dp_oc)
    dp_sc = np.array(dp_sc)
    white_refs = np.array(white_refs)
    
    # Calculate number of suns:
    num_suns = np.array(fluxes) / flux_1sun
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    rrs = (dp_oc-dp_sc)/(white_refs)
    Vs =  voc_rads +  (sci.k*298/sci.e)*np.log(rrs)
    Js = (1-num_suns)
    
    return num_suns, rrs, Vs, Js, fluxes


def whole_image_jv(path, flux_1sun, voc_rad, white_mean_scaled):
    # Lists to hold data
    dp_oc = []
    dp_sc = []
    white_refs = []
    fluxes = []
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
        
    filenames_voc = find_npy(f"{path}\\oc")    


    for filename_voc in filenames_voc:
        # Asuming there is a sc file for every oc file, and it's saved by the same name format:
        filename_isc = 'SC_'+'_'.join(filename_voc.split('_')[1:])
        flux = float(filename_voc.split('_')[1])
        exposure = float(filename_voc.split('_')[2])
        fluxes.append(flux)

        datapoint_voc =  np.mean(np.load(f"{path}\\oc\\{filename_voc}"))
        datapoint_isc =  np.mean(np.load(f"{path}\\sc\\{filename_isc}"))
        
        dp_oc.append(datapoint_voc/exposure)
        dp_sc.append(datapoint_isc/exposure)
        white_refs.append(white_mean_scaled*flux)

    dp_oc = np.array(dp_oc)
    dp_sc = np.array(dp_sc)
    white_refs = np.array(white_refs)
    
    # Calculate number of suns:
    num_suns = np.array(fluxes) / flux_1sun
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    rrs = (dp_oc-dp_sc)/(white_refs)
    Vs =  voc_rads +  (sci.k*298/sci.e)*np.log(rrs)
    Js = (1-num_suns)
    
    return num_suns, rrs, Vs, Js, fluxes


def array_jv(path, flux_1sun, voc_rad, white_mean_scaled):
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad/(sci.k*298))-1
        
    filenames_voc = find_npy(f"{path}\\oc") 

    for filename_voc in filenames_voc:
        # Asuming there is a sc file for every oc file, and it's saved by the same name format:
        filename_isc = 'SC_'+'_'.join(filename_voc.split('_')[1:])
        flux = float(filename_voc.split('_')[1])
        exposure = float(filename_voc.split('_')[2])

        im_voc =  np.load(f"{path}\\oc\\{filename_voc}")
        im_isc =  np.load(f"{path}\\sc\\{filename_isc}")
        
        white_ref = white_mean_scaled*flux 

        num_sun = flux/flux_1sun

    # Calculate number of suns:
    white_nonunif = np.load(f"{path}\\white.npy")
    num_suns = np.array(fluxes) / flux_1sun
    voc_rads = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_suns)+1)

    rrs = (dp_oc-dp_sc)/(white_refs)
    Vs =  voc_rads +  (sci.k*298/sci.e)*np.log(rrs)
    Js = (1-num_suns)
    
    return num_suns, rrs, Vs, Js, fluxes







