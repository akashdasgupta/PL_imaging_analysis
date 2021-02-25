from external_imports import *
from image_process import *

def single_pixel_jv(path, row, col, jsc, voc_rad,
                    nominal_vs, solar_scale_list, int_nonunif, 
                    white_imarr, white_voltage, white_scale_factor, ref):
    rrs = []
    num_suns = []
    Js = []
    Vs = []

    filenames = find_tif(path)

    for k in range(len(nominal_vs)):
        datapoint =  np.array(Image.open(path+'\\'+filenames[k]))[row,col] - ref[row,col]
        num_sun = solar_scale_list[k] * int_nonunif[row, col]
        if num_sun > 1:
            continue
        rr = datapoint * white_scale_factor / (white_imarr[row, col]*(ledf(nominal_vs[k])/ledf(white_voltage)))

        J = jsc * num_sun
        V = voc_rad + (sci.k*298/sci.e)*np.log(rr)

        rrs.append(rr)
        num_suns.append(num_sun)
        Js.append(J)
        Vs.append(V)
    return np.array(num_suns), np.array(rrs), np.array(Vs), np.array(Js)

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









