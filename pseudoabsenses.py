from geopy.distance import geodesic
from typing import Union, List

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from pyspatialml import Raster


class PseudoAbsences(object):
    def __init__(self, raster: Raster, xys: np.array) -> None:
        """
        raster -- растровые данные

        xys -- Двумерный массив точек присутствия (Y, X)!

        """
        self.raster = raster
        self.xys = xys

    def _generate(self, bg_occ_count: int, random_state: int, filter_point, filter_SRE) -> pd.DataFrame:
        """
        Функция генерации точек
        """
        np.random.seed(seed=random_state)

        # create np array to store randomly sampled data
        valid_samples = np.zeros((0, self.raster.count))
        valid_coordinates = np.zeros((0, 2))

        # loop until target number of samples is satisfied
        satisfied = False

        n = bg_occ_count
        while satisfied is False:

            # generate random row and column indices
            Xsample = np.random.choice(range(0, self.raster.width), n)
            Ysample = np.random.choice(range(0, self.raster.height), n)

            # create 2d numpy array with sample locations set to 1
            sample_raster = np.empty((self.raster.height, self.raster.width))
            sample_raster[:] = np.nan
            sample_raster[Ysample, Xsample] = 1

            # get indices of sample locations
            rows, cols = np.nonzero(np.isnan(sample_raster) == False)

            # convert row, col indices to coordinates
            xy = np.transpose(
                rasterio.transform.xy(self.raster.transform, rows, cols))

            if filter_point:
                xy = filter_point(xy)

            # sample at random point locations
            samples = self.raster.extract_xy(xy, return_array=True)

            # append only non-masked data to each row of X_random
            samples = samples.astype("float32").filled(np.nan)
            invalid_ind = np.isnan(samples).any(axis=1)
            samples = samples[~invalid_ind, :]
            xy = xy[~invalid_ind, :]
            
            if filter_SRE:
                invalid = filter_SRE(samples)
                samples = samples[~invalid, :]
                xy = xy[~invalid, :]

            valid_samples = np.append(valid_samples, samples, axis=0)            
            valid_coordinates = np.append(valid_coordinates, xy, axis=0)

            # check to see if target_nsamples has been reached
            if len(valid_samples) >= bg_occ_count:
                satisfied = True
            else:
                n = bg_occ_count - len(valid_samples)

        # ==========================================================================
        gdf = pd.DataFrame(valid_samples, columns=self.raster.names)
        gdf["geometry_X"] = valid_coordinates[:, 0]
        gdf["geometry_Y"] = valid_coordinates[:, 1]
        # ==========================================================================
        return gdf
    
    def generate_random(self, bg_occ_count: int = 10000, random_state: int = 1, excluding: bool = True) -> pd.DataFrame:
        """
        Cлучайный отбор из всех точек в исследуемой области («случайные»)

        bg_occ_count -- Количество точек псевдоотсутствия

        random_state -- integer to use within random.seed.   

        excluding -- исключая доступные точки присутствия 
        """
        if excluding:
            x, y = np.transpose(rasterio.transform.xy(self.raster.transform, (0, 1), (0, 1)))
            exclusion_buffer = geodesic((x[1], x[0]), (y[1], y[0])).km 
            gdf = self._generate(bg_occ_count, random_state, lambda xy: self._pount_filter(xy, exclusion_buffer), None)
        else: 
            gdf = self._generate(bg_occ_count, random_state, None, None)

        return gdf

    def generate_SRE(self, bg_occ_count: int = 10000, random_state: int = 1, args_signs: Union[List[str], None] = None) -> gpd.GeoDataFrame:
        """
        Cлучайный выбор точек из всех точек за пределами подходящей области, 
        оцененной прямолинейной поверхностной оболочкой из образца присутствия 
        (модель оболочки поверхностного диапазона, использующая только данные 
        только о присутствии)
        
        bg_occ_count -- Количество точек псевдоотсутствия
        
        random_state -- integer to use within random.seed.

        args_signs -- array названий признаков, участвующих в отборе точек 
        """
        gdf = pd.DataFrame(self.raster.extract_xy(self.xys, return_array=True), columns=self.raster.names)
        signs = []
        if args_signs is None:
            iterabl =  self.raster.names 
        else:
            iterabl =  args_signs
        for column in iterabl:
            signs.append((self.raster.names.index(column), np.percentile(gdf[column], 25), np.percentile(gdf[column], 75)))
            
        def surface_range_envelope(samples: np.array) -> np.array: 
            """
            surface_range_envelope
            """
            invalid = []
            for samp in samples:
                if all(True if percent[1] <= samp[percent[0]] <= percent[2] else False for percent in signs):
                    invalid.append(True)
                else: 
                    invalid.append(False)
            
            return np.array(invalid)

        return self._generate(bg_occ_count, random_state, None, surface_range_envelope)

    def generate_far(self,  bg_occ_count: int = 10000, random_state: int = 1, exclusion_buffer: float = 111) -> pd.DataFrame:
        """
        Cлучайный выбор любой точки, расположенной по меньшей мере на один градус
        широты или долготы от любой точки присутствия (метод "1°far") 

        Cлучайный выбор любой доступной точки, расположенной по меньшей мере в
        двух градусах от любой точки присутствия (метод "2°far")

        bg_occ_count -- Количество точек псевдоотсутствия

        random_state -- integer to use within random.seed. 

        exclusion_buffer -- Минимальное расстояние между данными о присутствии и данными псевдоотсутствия (км)

        Расстояние между двумя параллелями (фактически — дуга меридиана), которые отличаются на 1° по широте, — примерно 111 км.
        """

        return self._generate(bg_occ_count, random_state, lambda xy: self._pount_filter(xy, exclusion_buffer), None)

    def _pount_filter(self, xy: np.array, exclusion_buffer) -> np.array:
        """
        Исключение из сгенерированных точек, близких к точкам присутствия
        """
        result = np.zeros((0, 2))
        for xy_p in xy:
            flag = True
            for xys_p in self.xys:
                if geodesic((xy_p[1], xy_p[0]), (xys_p[1], xys_p[0])).km <= exclusion_buffer:
                    flag = False
                    break
            if flag:
                result = np.append(result, [xy_p], axis=0)
        return result
        