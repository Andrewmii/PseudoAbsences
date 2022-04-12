import os

import pandas as pd
from pyspatialml import Raster

from pseudoabsenses import PseudoAbsences


BASE_DIR_PROJECT = r"D:\Python\elapid-main\BioData" 
DATA_FOLDER = r"features-data"
OCC_FOLDER = "output-data"

class Project(object):
    bg_occ_count: int  = 10000
    static_path: str = ""


def load_pount(project: Project):
    raster = Raster([os.path.join(BASE_DIR_PROJECT, project.static_path, DATA_FOLDER, x) for x in os.listdir(os.path.join(BASE_DIR_PROJECT, project.static_path, DATA_FOLDER))])
    df_point = pd.read_csv(r"D:\Python\elapid-main\BioData\selected_occ.csv")
    output = raster.extract_xy(df_point.reindex(columns=["Unnamed: 0", "species", "Longitude", "Latitude"]).iloc[:,2:].to_numpy())
    return output
    

def generate_bg_points(project: Project):
    raster = Raster([os.path.join(BASE_DIR_PROJECT, project.static_path, DATA_FOLDER, x) for x in os.listdir(os.path.join(BASE_DIR_PROJECT, project.static_path, DATA_FOLDER))])
    # df = raster.to_pandas(resampling='nearest')
    # df.to_csv(os.path.join(BASE_DIR_PROJECT, project.static_path, OCC_FOLDER, '_bg_points_1.csv'), index=False) 
    df_point = pd.read_csv(r"D:\Python\elapid-main\BioData\selected_occ.csv")
    xys = df_point.reindex(columns=["Unnamed: 0", "species", "Longitude", "Latitude"]).iloc[:,2:].to_numpy()
    pseudo = PseudoAbsences(raster, xys)
    # pseudoabsence_points = pseudo.generate_far(bg_occ_count=60, random_state=1)
    pseudo.generate_far(bg_occ_count=5000, random_state=1).to_csv(os.path.join(BASE_DIR_PROJECT, project.static_path, OCC_FOLDER, 'bg_points_far.csv'), index=False) 
    pseudo.generate_random(bg_occ_count=5000, random_state=1).to_csv(os.path.join(BASE_DIR_PROJECT, project.static_path, OCC_FOLDER, 'bg_points_random.csv'), index=False) 
    pseudo.generate_SRE(bg_occ_count=5000, random_state=1).to_csv(os.path.join(BASE_DIR_PROJECT, project.static_path, OCC_FOLDER, 'bg_points_SRE.csv'), index=False) 
    # pseudo.generate_SRE(bg_occ_count=3000, random_state=1, args_signs=['bio_1', 'bio_12']).to_csv(os.path.join(BASE_DIR_PROJECT, project.static_path, OCC_FOLDER, 'bg_points_SRE.csv'), index=False) 
    # print(pseudoabsence_points)
    # pseudoabsence_points.to_csv(os.path.join(BASE_DIR_PROJECT, project.static_path, OCC_FOLDER, 'bg_points_3.csv'), index=False) 


if __name__ == "__main__":
    # load_pount(Project).to_csv(os.path.join(BASE_DIR_PROJECT, Project.static_path, OCC_FOLDER, 'points.csv'), index=False) 
    generate_bg_points(Project)