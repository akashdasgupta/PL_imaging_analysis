# PL IMAGING: ANALYSIS SCRIPT
*Analysis script for data collected using PL imaging setup, publication DOI: [WHEN AVALIABLE]*

*Developer: Akash Dasgupta*

*Scientific contributes: Suhas Mahesh, Pietro Caprioglio, Yen-Hung Lin, Karl-Augustin Zaininger, Robert D.J. Oliver, Philippe Holzhey, Suer Zhou, Melissa M. McCarthy, Joel A. Smith, Maximilian Frenzel1, M. Greyson Christoforo, James M. Ball, Bernard Wenger and Henry J. Snaith*

* If you use this code, please cite our paper in your publication: [DOI WHEN AVALIABLE]
* For data acquired using script: [LINK WHEN AVALIABLE]
* The analysis is expecting the raw data to be in a specific naming format, if you use your own scripts for acquisition make sure you emulate naming conventions of above

## Usage
* All functions stored in .py files, calibration data in the separate folder
* User may access these via jupyter notebooks
* Instructions for use are detailed in notebook markdown

* Notebooks provided are:
	* **Pre-process images.ipynb**: This takes very large raw tiff fata and saves it as in a cropped, averaged, npy array format. You will need to run this first before using the other notebooks, they are not designed to deal with raw data. 
	* **INTENSITY_DEPENDANT.ipyb**series_resistance.ipynb: Analysis of images taken at open and short circuit, intensity dependant
	* **VOLTAGE DEPENDANT.ipyb**: Analysis for images taken at constant illumination, sweeping voltage
	* **Imaging_Paper_figures_local.ipynb**: Used to create transport layer figures and cell figures for thr paper ([DOI WHEN AVALIABLE])
	* **series_resistance.ipynb**: DON'T USE THIS, IT'S NOT BEEN VALIDATED YET!!!
	
For queries please email: akash.dasgupta@physics.ox.ac.uk

