from external_imports import *
from image_process import *

def single_pixel_jv(pathvoc, pathisc, row, col, voc_rad,
                    nominal_vs, solar_scale_list, int_nonunif, 
                    white_imarr, white_voltage, white_scale_factor, cell_exposure_lists, 
                    ref_voc, ref_isc, ave_rad=5):
    PLQEs = []
    num_suns = []
    Vs = []

    filenames_voc = find_tif(pathvoc)    
    filenames_isc = find_tif(pathisc)
    
    for k in range(len(nominal_vs)):
        for candidate_filename in filenames_voc:
            if float(candidate_filename.split('_')[1].split('=')[1]) == nominal_vs[k]:
                filename_voc = candidate_filename
        for candidate_filename in filenames_isc:
            if float(candidate_filename.split('_')[1].split('=')[1]) == nominal_vs[k]:
                filename_isc = candidate_filename

        datapoint_voc =  np.mean(np.array(Image.open(pathvoc+'\\'+filename_voc))[row-ave_rad:row+ave_rad,col-ave_rad:col+ave_rad]- ref_voc[row-ave_rad:row+ave_rad,col-ave_rad:col+ave_rad])
        datapoint_isc =  np.mean(np.array(Image.open(pathisc+'\\'+filename_isc))[row-ave_rad:row+ave_rad,col-ave_rad:col+ave_rad]- ref_isc[row-ave_rad:row+ave_rad,col-ave_rad:col+ave_rad])
        
        datapoint = datapoint_voc-datapoint_isc
        num_sun = solar_scale_list[k] * np.mean(int_nonunif[row-ave_rad:row+ave_rad,col-ave_rad:col+ave_rad])
        PLQE = datapoint * (white_scale_factor/cell_exposure_lists[k]) / (np.mean(white_imarr[row-ave_rad:row+ave_rad,col-ave_rad:col+ave_rad])*(ledf(nominal_vs[k])/ledf(white_voltage)))
        V = voc_rad + (sci.k*298/sci.e)*np.log(PLQE)

        PLQEs.append(PLQE)
        num_suns.append(num_sun)
        Vs.append(V)

    return np.array(num_suns), np.array(PLQEs), np.array(Vs)

def array_jv(path, rmin, rmax, cmin, cmax, jsc, voc_rad,
                    nominal_vs, solar_scale_list, int_nonunif, 
                    white_imarr, white_voltage, white_scale_factor, ref, savepath):
    filenames = find_tif(path)
    
    # Make save dirs:
    if not os.isdir(savepath+'\\rr'):
        os.mkdir(savepath+'\\rr')
    if not os.isdir(savepath+'\\num_suns'):
        os.mkdir(savepath+'\\num_suns')
    if not os.isdir(savepath+'\\Js'):
        os.mkdir(savepath+'\\Js')
    if not os.isdir(savepath+'\\Vs'):
        os.mkdir(savepath+'\\Vs')

    for k in range(len(nominal_vs)):
        data_arr =  np.array(Image.open(path+'\\'+filenames[k]))[rmin:rmax,cmin:cmax] - ref[rmin:rmax,cmin:cmax]
        num_sun = solar_scale_list[k] * int_nonunif[rmin:rmax,cmin:cmax]
        if num_sun > 1:
            continue
        rr =  data_arr * white_scale_factor / (white_imarr[rmin:rmax,cmin:cmax]*(ledf(nominal_vs[k])/ledf(white_voltage)))

        J = jsc * num_sun
        V = voc_rad + (sci.k*298/sci.e)*np.log(rr)

        with open(savepath+f'\\rr\\rr_{k}', 'w') as file:
            np.save(file, rr)
        with open(savepath+f'\\num_suns\\numsuns_{k}', 'w') as file:
            np.save(file, num_sun)
        with open(savepath+f'\\J\\J_{k}', 'w') as file:
            np.save(file, J)
        with open(savepath+f'\\V\\V_{k}', 'w') as file:
            np.save(file, V)









