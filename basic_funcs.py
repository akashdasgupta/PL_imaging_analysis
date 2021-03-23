from external_imports import *
from PIL import Image

def path_process(path):
    db = {}  
    for root, dirs, _ in os.walk(path):
        for dir in dirs:
            db[dir] = []
            for _, subdirs, _ in os.walk(f"{root}\\{dir}"):
                for subdir in subdirs:
                    try:
                        int(subdir)
                        db[dir].append(int(subdir))
                    except ValueError:
                        pass
                break
        break
    return db

def find_tif(datapath):
    raw_paths = []
    for _, _, files in os.walk(datapath):
        for file in files:
            if file.endswith(".tif"): # Only interested in tifs
                raw_paths.append(file)
        break
    return raw_paths

# def image_for_cropping(db, basepath):
#     images = []
#     names = []
#     corosponding_paths = []
#     for key in db.keys():
#         for pix in db[key]:
#             image_names = find_tif(basepath+'\\'+key+'\\'+str(pix)+'\\camera')
#             nominal_vs = [i.split('_')[1].split('=')[1] for i in image_names]
#             max_filename = image_names[np.argmax(nominal_vs)]
#             im = np.array(Image.open(basepath+'\\'+key+'\\'+str(pix)+'\\camera\\'+max_filename))
#             im_scaled = im *(2**(8)/np.max(im)) # 8bit 
#             name = key+', Pixel '+str(pix)

#             images.append(im_scaled.astype(np.uint8))
#             names.append(name)
#             corosponding_paths.append(basepath+'\\'+key+'\\'+str(pix)+'\\camera')
#     return images, names, corosponding_paths




