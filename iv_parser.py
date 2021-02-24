from external_imports import *

def GetNominalV(path):
    nominal_v = []
    with open(path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            nominal_v.append(float(row[0]))
    return np.array(nominal_v)

def get_sm_data(path, illuminated_area):
    v = []
    i = []

    with open(path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            v.append(float(row[0]))
            i.append(float(row[1]))
    i = np.array(i)
    v = np.array(v)

    i /= illuminated_area
    return v,i

def nominalV2sunscale(nominal_vs, bandgap, lef_int_func):
    one_sun = j1sunf(bandgap)
    scale_photons = lef_int_func(nominal_vs)
    return scale_photons/one_sun


