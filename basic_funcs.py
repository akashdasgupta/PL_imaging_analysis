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

def get_nominal_V(path):
    nominal_v = []
    with open(path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            try:
                nominal_v.append(float(row[0]))
            except ValueError:
                pass
    return np.array(nominal_v)


