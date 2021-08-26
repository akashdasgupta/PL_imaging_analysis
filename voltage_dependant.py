from external_imports import *
from image_process import *
from basic_funcs import *
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

def vsweep_array_QFLS(path, voc_rad0):
    # For Voc_rad calculatuion:
    jsc_by_j0 = np.exp(sci.e*voc_rad0/(sci.k*298))-1
    
    for direction in ["vsweep_f", "vsweep_b", "vsweep"]:
        if os.path.isdir(f"{path}\\PLQE_{direction}"):

            filenames = find_npy(f"{path}\\PLQE_{direction}") 

            if not os.path.isdir(f"{path}\\QFLS_{direction}"):
                os.makedirs(f"{path}\\QFLS_{direction}")

            for filename in filenames:
                # Asuming there is a sc file for every oc file, and it's saved by the same name format:
                bias = float(filename.split('_')[0]) # should work for vsweep stuff
                flux = float(filename.split('_')[1])
                num_sun = float(filename.split('_')[1])

                PLQE =  np.load(f"{path}\\PLQE_{direction}\\{filename}")
                voc_rad = (sci.k*298/sci.e)*np.log((jsc_by_j0*num_sun)+1)

                QFLS =  voc_rad +  (sci.k*298/sci.e)*np.log(PLQE)
                # convention: [num suns]_[J]_[flux]_.npy
                savename = f"{bias}_{num_sun}_{flux}_"
                np.save(f"{path}\\QFLS_{direction}\\{savename}", QFLS)


def vsweep_col_eff(path, rawpath):
    for direction in ["vsweep_f", "vsweep_b", "vsweep"]:
        if not os.path.isdir(f"{path}\\{direction}_coleff"):
            os.makedirs(f"{path}\\{direction}_coleff")
        
        # Figure out if we have enough points for a OC image from this dataset
        Vs = []
        Is = []
        with open(f"{rawpath}\\{direction}\\source_meter.csv", 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                Vs.append(float(row[0]))
                Is.append(float(row[1]))
        
        sign0 = Is[0]/abs[Is[0]]
        oc_extrapolate = False
        for i in Is:
            if i/ abs(i) != sign0:
                oc_extrapolate = True
                V_above = Vs[i-1]
                I_above = Is[i-1]
                V_bellow = Vs[i]
                I_bellow = Is[i]
                break

        if oc_extrapolate:
            filenames = find_npy(f"{path}\\PLQE_{direction}") 
            for filename in filenames:
                if float(filename.split('_')[0].split('=')[1]) == np.around(V_above,3):
                    PLQE_oc_above = np.load(f"{path}\\PLQE_{direction}\\{filename}")
                elif float(filename.split('_')[0].split('=')[1]) == np.around(V_bellow, 3):
                    PLQE_oc_bellow = np.load(f"{path}\\PLQE_{direction}\\{filename}")
            
                m_ocarr = (PLQE_oc_bellow - PLQE_oc_above)/ (I_bellow-I_above)
                PLQE_oc = m_ocarr*I_bellow  # I=0, so oc image = c
        
        else:
            PLQE_oc = oc_image_maker(f"{path}\\PLQE_oc",f"{path}\\{direction}")
            if not PLQE_oc:
                raise FileNotFoundError("Couldn't find enough files for OC image")
        
        for filename in filenames:
            bias = float(filename.split('_')[0].split('=')[1])
            num_sun = float(filename.split('_')[1])
            flux = float(filename.split('_')[2])
            savename = f"{bias}_{num_sun}_{flux}_"

            PLQE_v =  np.load(f"{path}\\PLQE_{direction}\\{filename}")
            col_eff = (PLQE_oc-PLQE_v) / PLQE_oc
            np.save(f"{path}\\{direction}_coleff\\{savename}", col_eff)


def oc_image_maker(oc_path, vsweep_path):
    oc_filenames = find_npy(oc_path)
    vsweep_filenames = find_npy(vsweep_path)

    flux_vsweep = float(vsweep_filenames[0].split('_')[2])
    oc_fluxes = [float(i.split('_')[2]) for i in oc_filenames]
    oc_fluxes.sort()
    if max(oc_fluxes) < flux_vsweep or min(oc_fluxes) > flux_vsweep:
        return False
    
    sign0 = (oc_fluxes[0] - flux_vsweep) / abs(oc_fluxes[0] - flux_vsweep)
    for i, filename in enumerate(oc_filenames):
        flux = float(filename.split('_')[2])
        if (flux - flux_vsweep) / abs(flux - flux_vsweep) != sign0:
            oc_flux_above = flux
            oc_flux_bellow = float(oc_filenames[i-1].split('_')[2])
            oc_PLQE_above = oc_filenames[i]
            oc_PLQE_bellow = oc_filenames[i-1]
    
    m_arr = (oc_PLQE_bellow-oc_PLQE_above)/(oc_flux_bellow-oc_flux_above)
    c_arr = oc_flux_above - m_arr*oc_PLQE_above

    return (m_arr * flux_vsweep) + c_arr



