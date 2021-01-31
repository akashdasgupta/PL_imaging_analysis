from external_imports import *

def path_process(path):
    db = {}  

    for root, dirs, _ in os.walk(path):
        for dir in dirs:
            db[dir] = []
            for _, subdirs, _ in os.walk(root+"\\"+dir):
                for subdir in subdirs:
                    try:
                        float(subdir)
                        db[dir].append(subdir)
                    except ValueError:
                        pass
                break
        break
    return db


