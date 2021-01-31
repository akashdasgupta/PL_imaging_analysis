from external_imports import *

def GetNominalV(path):
    nominal_v = []
    corosponding_I = []
    with open(path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            nominal_v.append(float(row[0]))
            corosponding_I.append(float(row[1]))
    return nominal_v, corosponding_I

def averager(paths):
    temp = np.array(Image.open(paths[0]))
    rows, cols = temp.shape
    del temp

    buffer = np.zeros((rows, cols))
    for path in paths:
        im = np.array(Image.open(path), dtype=np.float32) # So there's enough headroom for the uint16
        buffer += im
        del im # unnescesary but just in case
    average_arr = buffer / len(paths)
    del buffer # unnescesary but just in case
    average_arr = average_arr.astype("uint16")
    return average_arr

def find_tif(datapath):
    raw_paths = []
    for _, _, files in os.walk(datapath):
        for file in files:
            if file.endswith(".tif"): # Only interested in tifs
                raw_paths.append(datapath + '/'+file)
    return raw_paths

def white_mask_maker(imarr):
    rows, cols = imarr.shape
    white_blurred = nd.gaussian_filter(imarr, 10)
    thresholded_arr = np.zeros((white_blurred.shape))
    half_max = np.max(white_blurred)/2.1
    for i in range(rows):
        for j in range(cols):
            if white_blurred[i,j] > half_max:
                thresholded_arr[i,j] = 0
            else:
                thresholded_arr[i,j] = 1
    rmin, cmin = (0,0)
    rmax,cmax = (rows,cols)
    for i in range(rows):
        if np.sum(thresholded_arr[i,:]) < cols:
            rmin = i
            break
    for j in range(cols):
        if np.sum(thresholded_arr[:,j]) < rows:
            cmin = j
            break
    for i in range(rows):
        if np.sum(thresholded_arr[rows-1-i,:]) < cols:
            rmax =rows-1-i
            break
    for j in range(cols):
        if np.sum(thresholded_arr[:,cols-1-j]) < rows:
            cmax = cols-1-j
            break
    
    return (rmin,rmax,cmin,cmax)

def calibrate(imarr, wimarr):
    return (imarr/wimarr )*np.mean(wimarr)




