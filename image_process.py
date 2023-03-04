# Parallel_processing
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

from external_imports import *
from basic_funcs import *


def calculate_spectral_correction(bandgap, filterODf):
    wavel_range =  np.arange(300, 1000, 1)
     # Get spectra for all involved components
    led_spec = led_specf(wavel_range)
    cell_spec = np.vectorize(BBf_cellf)(wavel_range, bandgap)
    cam_QE = camqef(wavel_range)
    lens_cal = lenscalf(wavel_range)
    filter_OD = filterODf(wavel_range)

    # Calculate relative factors
    white_factor =  integrate.simps((led_spec*cam_QE*lens_cal), wavel_range*1e-9) / integrate.simps(led_spec, wavel_range*1e-9)
    white_factor *= 0.99 # assume 99% reflectivity of ref 00
    cell_factor = integrate.simps((cell_spec*cam_QE*lens_cal*filter_OD), wavel_range*1e-9) / integrate.simps(cell_spec, wavel_range*1e-9)
    
    correction = white_factor/cell_factor
    return correction



def PLQEmap(datapath, filename, white_filepath, white_params_filepath, correction, flux_1sun):
    # Parse parameters
    bias = filename.split('_')[0]
    if bias.lower() != 'oc' and  bias.lower() != 'sc':
        bias = float(bias.split('=')[1])
    flux = float(filename.split('_')[1])
    num_sun = flux/flux_1sun

    # Load the white params: 
    with open(white_params_filepath,'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] == 'Exposure':
                white_exposure = float(row[1])
            if row[0] == 'Flux scale':
                white_flux = float(row[1])

    # load white
    white_mean = np.mean(np.load(white_filepath)) / white_exposure

    # Extract params from filename
    exposure = float(filename.split('_')[2])


    # Main calculation
    im_cell =  np.load(f"{datapath}/{filename}")/exposure
    PLQE = (im_cell/white_mean) * correction * (white_flux/flux)

    return bias, num_sun, flux, PLQE
