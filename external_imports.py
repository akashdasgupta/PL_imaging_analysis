import os
import csv
import numpy as np
from matplotlib import pyplot as plt
from scipy import constants as sci
from scipy.interpolate import interp1d as inter
import scipy.integrate as integrate
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects

#######################################################################################
# camera pixel to area at focal plane:

with open(r"./calibration_data/camera_pixel_area.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        pix = float(row[0])
        area = float(row[1])

pixels_per_cm = (pix/area)**0.5 # root of pixel area/ known area from cell pic
photodiode_diam_pix = pixels_per_cm*1.2
del pix
del area
#######################################################################################
# LED powermeter calibration:

nominal_v_cal = []
num_photons = []
with open(r"./calibration_data/ledcal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        nominal_v_cal.append(float(row[0]))
        num_photons.append(float(row[1]))
ledf = inter(nominal_v_cal,num_photons) # function
del nominal_v_cal # Free up memory
del num_photons
#######################################################################################
# LED spectrum (From Thorlabs):

led_spec_wavel = []
led_spec = []
with open(r"./calibration_data/ledwavel.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        led_spec_wavel.append(float(row[0]))
        led_spec.append(float(row[1]))
ledspecf = inter(led_spec_wavel, led_spec)
del led_spec_wavel
del led_spec
#######################################################################################
# Room temperature black body:

BBf = lambda wavel_bb: (100 * (560 / wavel_bb) ** 5 * (((np.exp((1.435 * 10 ** 7) 
                        / (298 * 560)) - 1) / (np.exp((1.435 * 10 ** 7) / 
                        (298 * wavel_bb)) - 1)))) 
def BBf_cellf(wavel, bandgap): 
    if (sci.h*sci.c/(wavel*1e-9*sci.e) > bandgap):
        return BBf(wavel) 
    else:
        return 0
#######################################################################################
# Camera Quantum efficiency: 

cam_qe_wavel = []
cam_qe = []
with open(r"./calibration_data/camqe.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        cam_qe_wavel.append(float(row[0]))
        cam_qe.append(float(row[1]))
camqef = inter(cam_qe_wavel, cam_qe)
del cam_qe_wavel
del cam_qe
#######################################################################################
# Filter response function:

filt_cal_wavel = []
filt_cal = []
with open(r"./calibration_data/filtcal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        filt_cal_wavel.append(float(row[0]))
        filt_cal.append(float(row[1]))
filtcalf = inter(filt_cal_wavel, filt_cal)
del filt_cal_wavel
del filt_cal

#######################################################################################
# Filter response function 2:

filt_cal_wavel = []
filt_cal = []
with open(r"./calibration_data/filtcal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        filt_cal_wavel.append(float(row[0]))
        filt_cal.append(float(row[1])**2)
filtcalf2 = inter(filt_cal_wavel, filt_cal)
del filt_cal_wavel
del filt_cal

#######################################################################################
# Lens calibration (alledgely):

lens_cal_wavel = []
lens_cal = []
with open(r"./calibration_data/lenscal.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        lens_cal_wavel.append(float(row[0]))
        lens_cal.append(float(row[1]))
lenscalf = inter(lens_cal_wavel, lens_cal)
del lens_cal_wavel
del lens_cal
#######################################################################################
# J_1sun (from AM1.5 spectrup)

bandgap_i = []
am_1_5 = []
with open(r"./calibration_data/AM_1_5_photon.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bandgap_i.append(float(row[0]))
        am_1_5.append(float(row[1]))
am15f = inter(bandgap_i, am_1_5)
del bandgap_i
del am_1_5
j1sunf = lambda bandgap : integrate.simps(am15f(np.arange(bandgap,4,0.001)),
                                          np.arange(bandgap,4,0.001))
#######################################################################################
# Voc, radiative limit:

bandgap_v = []
voc_rad = []
with open(r"./calibration_data/voc_rad_eg.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bandgap_v.append(float(row[0]))
        voc_rad.append(float(row[1]))
vocradf = inter(bandgap_v, voc_rad)
del bandgap_v
del voc_rad

#######################################################################################
# Area for 1 device (dark area on solar sim)
with open(r"./calibration_data/device_area.csv", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        device_area = float(row[1])

