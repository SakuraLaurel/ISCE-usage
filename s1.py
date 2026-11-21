import numpy as np
from lxml import etree
from os.path import join
from ast import literal_eval
from datetime import datetime
from functions import find_files, lreduce


class Algorithm(object):
    a = 6378137.0000  # 地球半长轴
    b = 6356752.3141  # 地球半短轴
    c = 299792458  # 光速，不是地球焦半径
    e2 = 1 - (b / a) ** 2  # 偏心率的平方
    ep2 = (a / b) ** 2 - 1  # 第二偏心率的平方
    f = 1 - b / a  # 扁率
    UTMScaleFactor = 0.9996  # UTM坐标系中央经线的比例因子
    n = (a - b) / (a + b)  # 我不知道这叫啥，但在坐标转换的时候要用

    @classmethod
    def polint(cls, xa, ya, x):
        """
        Neville插值法，详见 http://phyweb.physics.nus.edu.sg/~phywjs/CZ5101/NR-lec3.pptx
        params:
        xa: state_vectors对应的慢时间
        ya: state_vectors的某一个分量，如px，或者py、pz
        x: 要插值的时间
        """
        ns = np.argmin(np.abs(xa - x), axis=0)
        n, nx = xa.shape
        col_index = np.arange(nx)
        c, d, y = ya.copy(), ya.copy(), ya[ns, col_index]
        ns -= 1
        for m in range(1, n):
            for i in range(0, n - m):
                ho = xa[i] - x
                hp = xa[i + m] - x
                w = c[i + 1] - d[i]
                den = ho - hp
                # if np.any(abs(den) < abs(ho) * 1e-4):
                #     raise Exception('ERROR: subroutine polint, infinite slope!')
                den = w / den
                d[i] = hp * den
                c[i] = ho * den
            less = 2 * (ns + 1) < n - m
            y[less] += c[ns + 1, col_index][less]
            y[~less] += d[ns, col_index][~less]
            ns[~less] -= 1
        return y

    @classmethod
    def sod(cls, s: str, format: str):
        """
        日积秒
        """
        time = datetime.strptime(s, format)
        reference = time.replace(hour=0, minute=0, second=0, microsecond=0)
        return (time - reference).total_seconds()

    @classmethod
    def wgs2ecs(cls, lon, lat, height):
        """
        经纬度转地心直角坐标

        :param lon: 经度，角度制
        :param lat: 纬度，角度制
        :param height: 高程，米
        :return: [x, y, z] (ndarray, 单位 米)
        """
        lon = Algorithm.degree2rad(lon)
        lat = Algorithm.degree2rad(lat)
        N = Algorithm.a / np.sqrt(1 - Algorithm.e2 * np.sin(lat) ** 2)
        x = (N + height) * np.cos(lat) * np.cos(lon)
        y = (N + height) * np.cos(lat) * np.sin(lon)
        z = (N * (1 - Algorithm.e2) + height) * np.sin(lat)
        return np.array((x, y, z))

    @classmethod
    def ecs2wgs(cls, x, y, z):
        """
        地心直角坐标转经纬度

        :param x: x坐标值，米
        :param y: y坐标值，米
        :param z: z坐标值，米
        :return: [lon, lat, height] （经度，纬度，高程）（角度制，米）
        """
        p = np.sqrt(x**2 + y**2)
        theta = np.arctan(z * Algorithm.a / (p * Algorithm.b))
        lon = np.arctan2(y, x)
        lat = np.arctan(
            (z + Algorithm.e2 / (1 - Algorithm.e2) * Algorithm.b * np.sin(theta) ** 3)
            / (p - Algorithm.e2 * Algorithm.a * np.cos(theta) ** 3)
        )
        N = Algorithm.a / np.sqrt(1 - Algorithm.e2 * np.sin(lat) ** 2)
        height = p / np.cos(lat) - N
        lon, lat = Algorithm.rad2degree(lon), Algorithm.rad2degree(lat)
        return np.array((lon, lat, height))

    @classmethod
    def degree2rad(cls, degree):
        """
        角度制转弧度制

        :param degree: 角度值
        :return: 弧度值
        """
        return degree * np.pi / 180

    @classmethod
    def rad2degree(cls, rad):
        """
        弧度制转角度制

        :param rad: 弧度值
        :return: 角度值
        """
        return rad * 180 / np.pi

class Dem(object):
    def __init__(self) -> None:
        path = find_files("demLat_N38_N41_Lon_E114_E118.dem")[0]
        tree = etree.parse("%s.xml" % path)
        xp = tree.xpath("/demimage_name/component[@name='coordinate1']")[0]
        self.left = np.float32(xp.xpath("property[@name='startingvalue']/value")[0].text)
        self.right = np.float32(xp.xpath("property[@name='endingvalue']/value")[0].text)
        self.width = np.int32(xp.xpath("property[@name='size']/value")[0].text)
        yp = tree.xpath("/demimage_name/component[@name='coordinate2']")[0]
        self.top = np.float32(yp.xpath("property[@name='startingvalue']/value")[0].text)
        self.bottom = np.float32(yp.xpath("property[@name='endingvalue']/value")[0].text)
        self.length = np.int32(yp.xpath("property[@name='size']/value")[0].text)
        self.data = np.fromfile(path, dtype=np.int16).reshape((self.length, self.width)) - 10 # 粗略弥补没有EGM的不足
        self.dy = (self.top - self.bottom) / self.length
        self.dx = (self.right - self.left) / self.width


    def find(self, lon, lat):
        i = np.int32((self.top - lat) / self.dy)
        j = np.int32((lon - self.left) / self.dx)
        if any((i < 0) | (j < 0)):
            raise ValueError("超出DEM范围")
        return np.float32(self.data[i, j])

class Burst(object):
    def __init__(self, path, swath, burst, flat_phase=False):
        xml_path = join(path, "fine_coreg/IW%d.xml" % swath)
        xpath = (
            "/productmanager_name/component[@name='instance']/component[@name='bursts']/component[@name='burst%d']"
            % burst
        )
        tree = etree.parse(xml_path).xpath(xpath)[0]
        p = tree.xpath(
            "component[@name='orbit']/component[@name='state_vectors']/component"
        )
        time: list[datetime] = []
        self.state_vector_position = np.empty((len(p), 3), dtype=np.float32)
        self.state_vector_velocity = np.empty((len(p), 3), dtype=np.float32)
        for i, j in enumerate(sorted(p, key=lambda x: int(x.attrib["name"][11:]))):
            time.append(
                Algorithm.sod(
                    j.xpath("property[@name='time']/value")[0].text,
                    "%Y-%m-%d %H:%M:%S.%f",
                )
            )
            self.state_vector_position[i] = literal_eval(
                j.xpath("property[@name='position']/value")[0].text
            )
            self.state_vector_velocity[i] = literal_eval(
                j.xpath("property[@name='velocity']/value")[0].text
            )
        self.time_of_first_state_vector = time[0]
        self.state_vector_interval = time[1] - time[0]
        self.start_time = Algorithm.sod(
            tree.xpath("property[@name='sensingstart']/value")[0].text,
            "%Y-%m-%d %H:%M:%S.%f",
        )
        self.azimuth_line_time = float(
            tree.xpath("property[@name='azimuthtimeinterval']/value")[0].text
        )
        self.near_range_slc = float(
            tree.xpath("property[@name='startingrange']/value")[0].text
        )
        self.range_pixel_spacing = float(
            tree.xpath("property[@name='rangepixelsize']/value")[0].text
        )
        self.azimuth_lines = int(
            tree.xpath("property[@name='numberoflines']/value")[0].text
        )  # 1506
        self.range_samples = int(
            tree.xpath("property[@name='numberofsamples']/value")[0].text
        )  # 25441
        self.wavelength = float(
            tree.xpath("property[@name='radarwavelength']/value")[0].text
        )
        self.shape = (self.azimuth_lines, self.range_samples)
        self.slc_path = join(path, "fine_coreg/IW%d/burst_%02d.slc" % (swath, burst))
        self.rg_path = join(path, "fine_offsets/IW%d/range_%02d.off" % (swath, burst))
        self.lon_path = join(path, f"geom_reference/IW{swath}/lon_{burst:02d}.rdr")
        self.lat_path = join(path, f"geom_reference/IW{swath}/lat_{burst:02d}.rdr")
        self.flat_phase = flat_phase  # 是否要去除平地相位

    def _state(self, t, pos=True):
        # 计算卫星的位置、速度状态
        states = getattr(
            self, "state_vector_position" if pos else "state_vector_velocity"
        )
        n_states = 8
        initial = np.int16(
            np.around(
                (t - self.time_of_first_state_vector) / self.state_vector_interval
                - (n_states - 1) / 2
            )
        )
        up, down = len(states) - n_states, 0
        initial[initial < down] = down
        initial[initial > up] = up
        initial = np.expand_dims(initial, axis=0)
        index = np.repeat(initial, n_states, axis=0) + np.arange(n_states).reshape(
            (-1, 1)
        )
        xa = self.time_of_first_state_vector + index * self.state_vector_interval
        ya = states[index]
        x = Algorithm.polint(xa, ya[:, :, 0], t)
        y = Algorithm.polint(xa, ya[:, :, 1], t)
        z = Algorithm.polint(xa, ya[:, :, 2], t)
        return np.array((x, y, z))

    def load(self):
        # 在使用影像数据、地理坐标数据前必须调用本函数
        if hasattr(self, "data"):
            return
        self.data = np.fromfile(self.slc_path, np.complex64).reshape(self.shape)
        if not self.flat_phase:
            rg_offset = np.fromfile(self.rg_path, np.float32).reshape(self.shape)
            phase = np.exp(1j * 4 * np.pi * rg_offset * self.range_pixel_spacing / self.wavelength)
            self.data *= phase
        self.lon = np.fromfile(self.lon_path, np.float64).reshape(self.shape)
        self.lat = np.fromfile(self.lat_path, np.float64).reshape(self.shape)

    def unload(self):
        self.data, self.lon, self.lat = None, None, None

    def indirect(self, lon, lat, times=5):
        not_array = type(lon) != np.ndarray
        if not_array:
            lon = np.ones((1,1)) * lon
            lat = np.ones((1,1)) * lat
        h = Dem().find(lon, lat)
        R_T = Algorithm.wgs2ecs(lon, lat, h).reshape((3, -1))
        t = self.start_time * np.ones(R_T.shape[1])
        for _ in range(times):
            V_S = self._state(t, False)
            V_S_norm = np.linalg.norm(V_S, axis=0)
            R_ST = R_T - self._state(t, True)
            delta_r = np.sum(np.multiply(R_ST, V_S), axis=0) / V_S_norm
            t += delta_r / V_S_norm

        row = (t - self.start_time - self.azimuth_line_time / 2) / self.azimuth_line_time
        col = (
            np.linalg.norm(R_T - self._state(t, True), axis=0)
            - self.near_range_slc
            - self.range_pixel_spacing / 2
        ) / self.range_pixel_spacing
        row = np.int32(np.around(row.reshape(np.shape(lon))))
        col = np.int32(np.around(col.reshape(np.shape(lon))))
        if not_array:
            return row[0, 0], col[0, 0]
        else:
            return row, col

    def direct(self, row, col):
        self.load()
        row, col = np.int32(np.around(row)), np.int32(np.around(col))
        return self.lon[row, col], self.lat[row, col]

    def get(self, top, left, bottom, right):
        return self.data[bottom: top+1, left: right+1]

    def deburst(self, latter):
        _, lat1 = self.direct(-1, 0)
        _, lat2 = latter.direct(0, 0)
        _, lat3 = latter.direct(2, 0)
        nlat = np.int32(np.around((lat1 - lat2) / (lat3 - lat2)))
        for j in range(nlat-5, nlat+6):
            for k in range(-j-6, -j+4):
                if np.abs(self.direct(k, 0)[0] - latter.direct(j, 0)[0]) < 1e-5:
                    for i in ("data", "lon", "lat"):
                        setattr(self, i, np.concat([getattr(self, i)[:k], getattr(latter, i)[j:]], axis=0))
                    return self
        raise ValueError("找不到匹配的点")
            

class Swath(object):
    def __init__(self, path, swath):
        self.bs = [Burst(path, swath, i) for i in range(1, 10)] 
    
    def geolocate(self, top, left, bottom, right, inclination=11, ratio=np.cos(40*np.pi/180)):
        b1, row_min = 0, -1
        while row_min == -1:
            row, _ = self.bs[b1].indirect(right, bottom)
            if row > 0 and row < self.bs[b1].azimuth_lines:
                row_min = row
            else:
                b1 += 1
        b2, row_max = b1, -1
        while row_max == -1:
            row, _ = self.bs[b2].indirect(left, top)
            if row > 0 and row < self.bs[b2].azimuth_lines:
                row_max = row
            else:
                b2 += 1
        img = lreduce(lambda i,j:[j.deburst(i)], [self.bs[i] for i in range(b2, b1-1,-1)])[0]
        row_max, _ = img.indirect(left, top)
        _, col_min = img.indirect(left, bottom)
        _, col_max = img.indirect(right, top)
        theta = inclination * np.pi / 180
        coe = ratio*np.sin(theta)/(ratio*np.sin(theta)+np.cos(theta))*(ratio-np.tan(theta))/(ratio-ratio*np.tan(theta)**2)
        rmax = np.int32(row_max - coe *(row_max - row_min)) 
        rmin = np.int32(row_min + coe *(row_max - row_min)) 
        coe = np.sin(theta)/(np.sin(theta)+ratio*np.cos(theta))*(1-ratio*np.tan(theta))/(1-np.tan(theta)**2)
        cmax = np.int32(col_max - coe *(col_max - col_min)) 
        cmin = np.int32(col_min + coe *(col_max - col_min)) 
        content = [[0 for _ in range(5)], [0 for _ in range(5)]]
        content[0][0], content[1][0] = img.direct(rmax, cmax)
        content[0][1], content[1][1] = img.direct(rmax, cmin)
        content[0][2], content[1][2] = img.direct(rmin, cmin)
        content[0][3], content[1][3] = img.direct(rmin, cmax)
        content[0][4], content[1][4] = content[0][0], content[1][0]
        return img.get(rmax, cmin, rmin, cmax+1), content