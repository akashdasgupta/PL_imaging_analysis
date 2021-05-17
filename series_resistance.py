from external_imports import *
from image_process import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

def series_resistance(path):
    biased_filename = find_npy(f"{path}\\sr")[0]
    sc_filenames = find_npy(f"{path}\\sc")

    biased_bias = float(biased_filename.split('_')[0].split('=')[1])
    biased_j = sci.e * float(biased_filename.split('_')[1])

    biased_imarr = np.load(f"{path}\\sr\\{biased_filename}")

    series_resistance = np.zeros(biased_imarr.shape)
    def single_pix(i,jmax): 
        for j in range(jmax):
            unbiased_fluxes = []
            unbiased_intensities= []
            for sc_filename in sc_filenames:
                unbiased_intensities.append(np.load(f"{path}\\sc\\{sc_filename}")[i,j])
                unbiased_fluxes.append(sci.e * float(sc_filename.split('_')[1]))
            f = inter(unbiased_intensities,unbiased_fluxes, bounds_error=False)
            
            biased_intensity = biased_imarr[i,j]
            unbiased_j = f(biased_intensity)

            series_resistance_ij = abs(biased_bias)/abs(unbiased_j-biased_j)
            series_resistance[i,j] = series_resistance_ij
            
    Parallel(n_jobs=num_cores, backend='threading')(delayed(single_pix)(row, series_resistance.shape[1]) for row in range(series_resistance.shape[0]))

    return series_resistance