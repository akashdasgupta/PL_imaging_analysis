import os
import csv
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image
from joblib import Parallel, delayed
import multiprocessing
from scipy import integrate
from numpy.lib.function_base import delete
from scipy import signal as sg
from scipy import ndimage as nd
from scipy.interpolate import interp1d as inter
from skimage.feature import peak_local_max as findpeaks
from interactivecrop.interactivecrop import main as crop
from copy import deepcopy
from scipy import constants as sci

nominal_v_cal = []
num_photons = []
with open("ledcal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        nominal_v_cal.append(float(row[0]))
        num_photons.append(float(row[1]))

ledf = inter(nominal_v_cal,num_photons)
del nominal_v_cal
del num_photons

led_spec_wavel = []
led_spec = []
with open("ledwavel.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        led_spec_wavel.append(float(row[0]))
        led_spec.append(float(row[1]))

ledspecf = inter(led_spec_wavel, led_spec)
del led_spec_wavel
del led_spec

BBf = lambda wavel_bb: (100 * (560 / wavel_bb) ** 5 * (((np.exp((1.435 * 10 ** 7) / (298 * 560)) - 1) / 
                            (np.exp((1.435 * 10 ** 7) / (298 * wavel_bb)) - 1)))) 

cam_qe_wavel = []
cam_qe = []
with open("camqe.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        cam_qe_wavel.append(float(row[0]))
        cam_qe.append(float(row[1]))

camqef = inter(cam_qe_wavel, cam_qe)
del cam_qe_wavel
del cam_qe

filt_cal_wavel = []
filt_cal = []
with open("filtcal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        filt_cal_wavel.append(float(row[0]))
        filt_cal.append(float(row[1]))

filtcalf = inter(filt_cal_wavel, filt_cal)
del filt_cal_wavel
del filt_cal


lens_cal_wavel = []
lens_cal = []
with open("lenscal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        lens_cal_wavel.append(float(row[0]))
        lens_cal.append(float(row[1]))

lenscalf = inter(lens_cal_wavel, lens_cal)
del lens_cal_wavel
del lens_cal