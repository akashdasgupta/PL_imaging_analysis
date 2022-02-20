from PIL import Image
from scipy import ndimage as nd # Gausian filter
from interactivecrop.interactivecrop import main as crop # Cropping window
from copy import deepcopy # So we don't screw up raw arrays that we act on
from scipy.ndimage import gaussian_filter

# Parallel_processing
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

from external_imports import *
from basic_funcs import *

def averager(basepath, filenames):
    '''
    Opens and returns average of tiff files
    '''
    temp = np.array(Image.open(basepath+'/'+filenames[0]))
    rows, cols = temp.shape
    del temp

    buffer = np.zeros((rows, cols))
    for filename in filenames:
        im = np.array(Image.open(basepath+'/'+filename), dtype=np.float32) # So there's enough headroom for the uint16
        buffer += im
        del im # unnescesary but just in case
    average_arr = buffer / len(filenames)
    del buffer # unnescesary but just in case
    return average_arr


def white_nonunif(imarr, center_row, center_col, diam):
    """
    Returns array corosponding to per-pixel scale factor for photodiode measured flux
    * Takes mean of area equal to photodiode area, and sets this mean as the photodiode measurement
    * Normalises then for nonuniformaty
    """
    pix_in_powermeter = []
    for i in range(imarr.shape[0]):
        for j in range(imarr.shape[1]):
            if( (i-center_row)**2 + (j-center_col)**2) <= (diam/2)**2:
                pix_in_powermeter.append(imarr[i,j])
    cal_arr = deepcopy(imarr)
    return cal_arr/ np.mean(pix_in_powermeter)

def beam_correct(imarr, white_imarr_norm, ref_imarr, ledv, rmin, rmax, cmin, cmax):
    """
    Corrects for effect of beam nonuniformity in the raw counts image
    """ 
    cropped_white_norm = white_imarr_norm[rmin:rmax, cmin:cmax]
    cropped_raw_image = imarr[rmin:rmax, cmin:cmax]
    cropped_ref_imarr = ref_imarr[rmin:rmax, cmin:cmax]

    overall_photon_flux = ledf(ledv)
    photon_flux_on_cell = np.mean(cropped_white_norm)*overall_photon_flux

    return photon_flux_on_cell, ((cropped_raw_image-cropped_ref_imarr)/cropped_white_norm)*np.mean(cropped_white_norm)

def white_over_cell_correction(led_specf, cell_specf, bandgap, camQEf, lenscalf, white_reflectivity, filterODf):
    wavel_range =  np.arange(300, 1000, 1)

    led_spec = led_specf(wavel_range)
    cell_spec = cell_specf(wavel_range, bandgap)

    cam_QE = camQEf(wavel_range)
    lens_cal = lenscalf(wavel_range)
    filter_OD = filterODf(wavel_range)

    white_factor =  integrate.simps((led_spec*cam_QE*lens_cal), wavel_range*1e-9) / integrate.simps(led_spec, wavel_range*1e-9)
    white_factor *= white_reflectivity 

    cell_factor = integrate.simps((cell_spec*cam_QE*lens_cal*filter_OD), wavel_range*1e-9) / integrate.simps(cell_spec, wavel_range*1e-9)

    return white_factor/cell_factor


def PLQEmap(datapath, filename, whitefilepath, whiteparamsfilepath, bandgap, filterODf, flux_1sun=None, one_sun_flux_pull_path=None):
    bias = filename.split('_')[0]
    if bias.lower() != 'oc' and  bias.lower() != 'sc':
        bias = float(bias.split('=')[1])
    flux = float(filename.split('_')[1])

    # Flux at 1sun
    if not flux_1sun :
        flux_1sun = j1sunf(bandgap)
    if one_sun_flux_pull_path:
        name = find_npy(one_sun_flux_pull_path)
        flux_1sun = float(name[0].split('_')[1])

    print(flux_1sun)
    num_sun = flux/flux_1sun

           
    # Load the white params: 
    white_exposure = None
    white_flux = None

    with open(whiteparamsfilepath,'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] == 'Exposure':
                white_exposure = float(row[1])
            if row[0] == 'Flux scale':
                white_flux = float(row[1])

    if white_exposure == None or white_flux == None:
        raise ValueError("Couldn't find white params!")

    # calculate correction
    correction = white_over_cell_correction(ledspecf, np.vectorize(BBf_cellf), 
                                                bandgap, camqef, lenscalf, 0.99, filterODf)

    # load white
    white_mean = np.mean(np.load(whitefilepath)) / white_exposure

    # Extract params from filename
    exposure = float(filename.split('_')[2])


    # Main calculation
    im_cell =  np.load(f"{datapath}/{filename}")/exposure
    PLQE = (im_cell/white_mean) * correction * (white_flux/flux)

    return bias, num_sun, flux, PLQE

def save_PLQE(datapath, savepath, whitefilepath, whiteparamsfilepath, bandgap, filterODf, flux_1sun=None, savename='untitled', one_sun_flux_pull_path=None):
    if not os.path.isdir(f"{savepath}/PLQE_{savename}"):
        os.makedirs(f"{savepath}/PLQE_{savename}")
    filenames = find_npy(datapath)
    def process_one_file(filename):
        bias, num_sun, flux, PLQE =  PLQEmap(datapath,filename, whitefilepath, whiteparamsfilepath, bandgap, filterODf, flux_1sun=flux_1sun,one_sun_flux_pull_path=one_sun_flux_pull_path)
        savefilename = f"{bias}_{num_sun}_{flux}_"
        np.save(f"{savepath}/PLQE_{savename}/{savefilename}", PLQE)
    
    Parallel(n_jobs=num_cores)(delayed(process_one_file)(filename) for filename in filenames)

ddef external_PLQE_averager(datapath, savepath):

    path_db = path_process(datapath)
    for key in path_db.keys():

        for pix in path_db[key]:
            for _,biases,_ in os.walk(f"{datapath}/{key}/{pix}"):
                break
            break
        
        for bias in biases:
            if bias.lower() == 'oc' or bias.lower() == 'sc':
                 # We assume the file holds the 1 sun value
                with open(f"{datapath}/{key}/PLQE_EXT/{bias}/plqe_ext.txt", 'r') as file:
                    reader = csv.reader(file)
                    for row in reader:
                        ext_PLQE_1sun = float(row[0])
                        break
            pix_ave_plqe = {}
            for pix in path_db[key]:
                filenames = find_npy(f"{savepath}/{key}/{pix}/PLQE_{bias}")
                numsuns = [float(i.split('_')[1]) for i in filenames]
                means = [np.mean(np.load(f"{savepath}/{key}/{pix}/PLQE_{bias}/{i}")) for i in filenames]

                means = [i for _,i in sorted(zip(numsuns, means))]
                numsuns.sort()
                
                if len(numsuns) != 1:
                    pix_ave_plqe[pix] = inter(numsuns, means, fill_value="extrapolate")(1)
                else:
                    pix_ave_plqe[pix] = means[0]

            correction = ext_PLQE_1sun / np.mean([pix_ave_plqe[i] for i in  pix_ave_plqe.keys()])
            print(pix_ave_plqe, ext_PLQE_1sun)



            for pix in path_db[key]:
                filenames = find_npy(f"{savepath}/{key}/{pix}/PLQE_{bias}")
                for filename in filenames:
                    arr = np.load(f"{savepath}/{key}/{pix}/PLQE_{bias}/{filename}")
                    corrected_arr = arr * correction
                    np.save(f"{savepath}/{key}/{pix}/PLQE_{bias}/{filename}", corrected_arr)



            else:
                # We assume the file holds the oc value
                pass # I'll do this later


def fftbin(imagepath):
    cart_im = np.load(imagepath)
            
    r_binwidth = 1.5
    r_max = np.sqrt(np.sum([i**2 for i in cart_im.shape]))/2
    
    r_binedges = np.arange(0, r_max, r_binwidth)
    
    
    theta_summed = np.zeros(r_binedges.shape)
    theta_num  = np.zeros(r_binedges.shape)
    np.save('temp1', theta_summed)
    np.save('temp2', theta_num)
    
    def process_row(k):
        theta_summed = np.load('temp1.npy', mmap_mode='r+')
        theta_num = np.load('temp2.npy', mmap_mode='r+')
        for i in range(cart_im.shape[0]):
            for j in range(cart_im.shape[1]):
                x = i - cart_im.shape[0]/2
                y = cart_im.shape[1]/2 - j
                r = np.sqrt(x**2 + y**2)

                r_binedge = r_binedges[k]
                if r_binedge <= r < r_binedge + r_binwidth:
                    theta_summed[k] += cart_im[i,j]
                    theta_num += 1
    
    
    Parallel(n_jobs=num_cores, verbose = 0)(delayed(process_row)(k) for k in range(len(r_binedges)))
    
    theta_summed_f = np.load('temp1.npy')
    theta_num_f = np.load('temp2.npy')
    os.remove('temp1.npy')
    os.remove('temp2.npy')
    
    
    return theta_summed_f / theta_num_f
            



white_buffer=[None, None, None, None]
def callback_return_size(image_name, shape):
    colstart, rowstart, width, height = shape.size
    rowend = rowstart+height
    colend = colstart+width
    
    white_buffer[0] = rowstart
    white_buffer[1] = rowend
    white_buffer[2] = colstart
    white_buffer[3] = colend


cell_borders={}
def callback_return_pix_size(image_name, shape):
    colstart, rowstart, width, height = shape.size
    rowend = rowstart+height
    colend = colstart+width
    
    pix_buffer=[]
    pix_buffer.append(rowstart)
    pix_buffer.append(rowend)
    pix_buffer.append(colstart)
    pix_buffer.append(colend)

    cell_borders[image_name] = pix_buffer

