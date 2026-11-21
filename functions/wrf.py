from functions.trans import rad2degree, cos_vectors
from functions.image import wgs2ecs, ecs2wgs
import xarray as xr
import numpy as np


class Wrf(object):
    """
    处理wrf运行结果
    """

    def __init__(self, path):
        data = xr.open_dataset(path)
        self.__set_params(data)
        self.__set_extent(data)
        self.__cal_N()
        self.__cal_i2h()

    def __set_params(self, data):
        self.p = np.array((data['P'][0] + data['PB'][0]) / 100)  # hPa
        self.h = np.array((data['PH'][0] + data['PHB'][0]) / 9.81)  # m
        self.vapor_pressure = np.array(
            self.p / (1 + 0.622 / data['QVAPOR'][0]))  # hPa
        hur = data['QVAPOR'][0]/(data['QVAPOR'][0]+1)
        ta = data['T'][0] + 300  # Kelvin
        self.t = np.array(ta*(self.p/1000) **
                          (287/1004*(1+0.61*hur)/(1+0.85*hur)))

    def __set_extent(self, data):
        def helper(center, resolution, n_edges):
            """
            返回的是边界上的起点和终点，也就是说比如中心坐标为(1,1)，分辨率为1，那么左右边界为0.5和1.5。
            但是n_edges-1是网格数，也就是说如果只有1个数，那么n_edges=2, 返回的n_edges-1为1。
            """
            radius = 6370e3  # wrf内部实现中，固定了半径是6370km，无论纬度是多少。
            resolution = rad2degree(resolution / radius)
            start = center - (n_edges / 2 - 0.5) * resolution
            end = start + (n_edges-1)*resolution
            return start, end, resolution, n_edges - 1

        self.left, self.right, self.dlon, self.nlon = helper(
            data.attrs['CEN_LON'], data.attrs['DX'], data.attrs['WEST-EAST_GRID_DIMENSION'])
        self.bottom, self.top, self.dlat, self.nlat = helper(
            data.attrs['CEN_LAT'], data.attrs['DY'], data.attrs['SOUTH-NORTH_GRID_DIMENSION'])

    def __cal_N(self):
        dry_pressure = self.p - self.vapor_pressure
        k1, k2, k3 = 77.6, 71.6, 3.75e5
        N_dry = k1 * dry_pressure / self.t
        N_wet = k2 * self.vapor_pressure / \
            self.t + k3 * self.vapor_pressure / self.t**2
        self.N = N_dry + N_wet
        self.offset = 1e-6*(k1*287.05*np.mean(self.p[-1])/9.81)

    def __cal_i2h(self):
        threshold_delta_pressure = 900
        split_index = np.argmax(
            (self.h[1:]-self.h[:-1]) > threshold_delta_pressure, axis=0)
        rows, cols = np.arange(self.h.shape[1]), np.arange(self.h.shape[2])
        cols, rows = np.meshgrid(cols, rows)
        h_left = self.h[split_index-1, rows, cols]
        h_right = self.h[split_index, rows, cols]
        rate = (self.h[1]-self.h[0]) / (h_right-h_left)
        log_rate = np.log(rate)
        coe1 = (h_left - self.h[0])/(1-rate)
        coe2 = h_left - coe1

        def func1(i, j, k):
            index = i*(-log_rate[j, k])/(split_index[j, k]-1)+log_rate[j, k]
            return coe1[j, k] * np.e**index+coe2[j, k]

        max_index = len(self.h) - 1
        coe3 = (self.h[max_index]-h_right)/(max_index-split_index)

        def func2(i, j, k):
            return coe3[j, k] * (i-split_index[j, k])+h_right[j, k]

        def i2h(i, j, k):
            i, j, k = np.array(i), np.array(j), np.array(k)
            res = np.empty(np.shape(i), dtype=np.float32)
            interval = i < split_index[j, k]
            res[interval] = func1(i[interval], j[interval], k[interval])
            res[~interval] = func2(i[~interval], j[~interval], k[~interval])
            return res
        self.i2h = i2h

        def h2i(h, j, k):
            res = np.empty(h.shape, dtype=np.float32)
            jj, kk = j, k
            interval = h < h_right[j, k]
            j, k = jj[interval], kk[interval]
            res[interval] = (np.log((h[interval] - coe2[j, k]) / coe1[j, k]) -
                             log_rate[j, k])*(split_index[j, k]-1)/-log_rate[j, k]
            j, k = jj[~interval], kk[~interval]
            res[~interval] = (h[~interval]-h_right[j, k]) / \
                coe3[j, k]+split_index[j, k]
            return res
        self.h2i = h2i

    def i2N(self, i, lat_is, lon_is):
        i -= 0.54
        i[i < 0] = 0
        i[i > len(self.N) - 1] = len(self.N) - 1
        i_down = np.int32(np.floor(i))
        i_up = np.int32(np.ceil(i))
        i_delta = i - i_down
        return self.N[i_down, lat_is, lon_is] * (1 - i_delta) + self.N[i_up, lat_is, lon_is] * i_delta

    def cal_delays(self, sp, lon, lat, h):
        if np.ndim(sp) == 1:  # 如果是一组数据
            res_shape = (1, )
            sp, lon, lat, h = sp.reshape((3, 1)), np.reshape(
                lon, (1,)), np.reshape(lat, (1,)), np.reshape(h, (1,))
        else:  # 如果是多组数据
            res_shape = sp.shape[1:]
            sp, lon, lat, h = sp.reshape(
                (3, -1)), lon.ravel(), lat.ravel(), h.ravel()
        
        tp = wgs2ecs(lon, lat, h)
        cos_zenith = cos_vectors(tp, sp - tp)
        dp = sp - tp
        coe = dp / (cos_zenith * np.linalg.norm(dp, axis=0))
        heights = self.get_heights(lon, lat)
        points = np.repeat(np.expand_dims(heights, 1), 3, 1) * coe.reshape((1,3, -1)) + tp.reshape((1,3, -1))
        
        lon, lat, h = ecs2wgs(points[:,0], points[:,1], points[:,2])
        lon_is = np.int32((lon - self.left) / self.dlon)[:-1]
        lat_is = np.int32((lat - self.bottom) / self.dlat)[:-1]
        indexes = self.h2i((heights[1:]+heights[:-1])/2, lat_is, lon_is)
        
        delays = 1e-6*np.sum(np.multiply(heights[1:]-heights[:-1], self.i2N(indexes, lat_is, lon_is)), axis=0)
        return delays.reshape(res_shape), cos_zenith.reshape(res_shape), np.ones(res_shape) * self.offset

    def get_heights(self, lon, lat):
        lon_is = np.int32((lon - self.left) / self.dlon)
        lat_is = np.int32((lat - self.bottom) / self.dlat)
        if np.any(lon_is < 0) or np.any(lon_is >= self.nlon) or np.any(lat_is < 0) or np.any(lat_is >= self.nlat):
            raise ValueError("超出WRF边界范围")
        return self.h[:, lat_is, lon_is]
    
    def get_params(self,param, lon, lat):
        lon_is = np.int32((lon - self.left) / self.dlon)
        lat_is = np.int32((lat - self.bottom) / self.dlat)
        if np.any(lon_is < 0) or np.any(lon_is >= self.nlon) or np.any(lat_is < 0) or np.any(lat_is >= self.nlat):
            raise ValueError("超出WRF边界范围")
        return param[:, lat_is, lon_is]