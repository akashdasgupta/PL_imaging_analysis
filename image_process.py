from PIL import Image
from scipy import ndimage as nd # Gausian filter
from interactivecrop.interactivecrop import main as crop # Cropping window
from copy import deepcopy # So we don't screw up raw arrays that we act on

from external_imports import *
from basic_funcs import *

def averager(basepath, filenames):
    temp = np.array(Image.open(basepath+'\\'+filenames[0]))
    rows, cols = temp.shape
    del temp

    buffer = np.zeros((rows, cols))
    for filename in filenames:
        im = np.array(Image.open(basepath+'\\'+filename), dtype=np.float32) # So there's enough headroom for the uint16
        buffer += im
        del im # unnescesary but just in case
    average_arr = buffer / len(basepath+'\\'+filename)
    del buffer # unnescesary but just in case
    return average_arr

def white_nonunif(imarr, center_row, center_col, diam):
    pix_in_powermeter = []
    for i in range(imarr.shape[0]):
        for j in range(imarr.shape[1]):
            if( (i-center_row)**2 + (j-center_col)**2) <= (diam/2)**2:
                pix_in_powermeter.append(imarr[i,j])
    cal_arr = deepcopy(imarr)
    return cal_arr/ np.mean(pix_in_powermeter)

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



'''
OLD FUNCTIONS: MAY USE SOME DAY LATER

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





'''