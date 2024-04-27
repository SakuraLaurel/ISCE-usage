from ast import literal_eval
from lxml import etree
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np


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


class Img(object):
    def __init__(self, xmlPath, slcPath):
        self._analyzeXML(xmlPath)
        self.data = np.fromfile(slcPath, np.complex64).reshape(
            (self.azimuth_lines, self.range_samples)
        )

    def geolocate(self, lon, lat, h, times=5):
        R_T = Algorithm.wgs2ecs(lon, lat, h).reshape((3, -1))
        t = self.start_time * np.ones(R_T.shape[1])
        for _ in range(times):
            V_S = self._state(t, False)
            V_S_norm = np.linalg.norm(V_S, axis=0)
            R_ST = R_T - self._state(t, True)
            delta_r = np.sum(np.multiply(R_ST, V_S), axis=0) / V_S_norm
            t += delta_r / V_S_norm

        row = (t - self.start_time) / self.azimuth_line_time
        col = (
            np.linalg.norm(R_T - self._state(t, True), axis=0) - self.near_range_slc
        ) / self.range_pixel_spacing
        row, col = row.reshape(np.shape(lon)), col.reshape(np.shape(lon))
        return row, col

    def find(self, y, x, threshold=3):
        slices = (
            (y < 0)
            | (y > self.azimuth_lines - 1)
            | (x < 0)
            | (x > self.range_samples - 1)
        )
        y[slices] = 0
        x[slices] = 0
        row, col = np.uint32(y), np.uint32(x)
        rr, cr = y - row, x - col
        v1 = (
            self.data[row, col],
            self.data[row, col + 1],
            self.data[row + 1, col],
            self.data[row + 1, col + 1],
        )
        v2 = ((1 - rr) * (1 - cr), (1 - rr) * cr, rr * (1 - cr), rr * cr)
        intensity = np.sum(np.multiply(np.abs(v1) ** 2, v2), axis=0, dtype=np.float32)
        intensity[slices] = np.nan
        max_v = np.nanpercentile(intensity, 100 - threshold)
        min_v = np.nanpercentile(intensity, threshold)
        intensity[intensity > max_v] = max_v
        intensity[intensity < min_v] = min_v
        return intensity

    def _state(self, t, pos=True):
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


class Stripmap(Img):
    def _analyzeXML(self, path):
        tree = etree.parse(path).xpath("/productmanager_name/component[@name='instance']")[0]
        p = tree.xpath("component[@name='orbit']/component[@name='state_vectors']/component")
        time: list[datetime] = []
        self.state_vector_position = np.empty((len(p), 3), dtype=np.float32)
        self.state_vector_velocity = np.empty((len(p), 3), dtype=np.float32)
        for i, j in enumerate(sorted(p, key=lambda x: int(x.attrib["name"][11:]))):
            time.append(Algorithm.sod(j.xpath("property[@name='time']/value")[0].text, "%Y-%m-%d %H:%M:%S"))
            self.state_vector_position[i] = literal_eval(j.xpath("property[@name='position']/value")[0].text)
            self.state_vector_velocity[i] = literal_eval(j.xpath("property[@name='velocity']/value")[0].text)
        self.time_of_first_state_vector = time[0]
        self.state_vector_interval = time[1] - time[0]
        self.start_time = Algorithm.sod(tree.xpath("property[@name='sensing_start']/value")[0].text, "%Y-%m-%dT%H:%M:%S.%f")
        self.azimuth_line_time = 1 / float(tree.xpath("component[@name='instrument']/property[@name='prf']/value")[0].text)
        self.near_range_slc = float(tree.xpath("property[@name='starting_range']/value")[0].text)
        self.range_pixel_spacing = float(tree.xpath("component[@name='instrument']/property[@name='range_pixel_size']/value")[0].text)
        self.azimuth_lines = int(tree.xpath("property[@name='number_of_lines']/value")[0].text)
        self.range_samples = int(tree.xpath("property[@name='number_of_samples']/value")[0].text)

class Tops(Img):
    def __init__(self, xmlPath, slcPath, burst):
        self.burst: int = burst
        super().__init__(xmlPath, slcPath)

    def _analyzeXML(self, path):
        tree = etree.parse(path).xpath("/productmanager_name/component[@name='instance']/component[@name='bursts']/component[@name='burst%d']" % self.burst)[0]
        p = tree.xpath("component[@name='orbit']/component[@name='state_vectors']/component")
        time: list[datetime] = []
        self.state_vector_position = np.empty((len(p), 3), dtype=np.float32)
        self.state_vector_velocity = np.empty((len(p), 3), dtype=np.float32)
        for i, j in enumerate(sorted(p, key=lambda x: int(x.attrib["name"][11:]))):
            time.append(Algorithm.sod(j.xpath("property[@name='time']/value")[0].text, "%Y-%m-%d %H:%M:%S.%f"))
            self.state_vector_position[i] = literal_eval(j.xpath("property[@name='position']/value")[0].text)
            self.state_vector_velocity[i] = literal_eval(j.xpath("property[@name='velocity']/value")[0].text)
        self.time_of_first_state_vector = time[0]
        self.state_vector_interval = time[1] - time[0]
        self.start_time = Algorithm.sod(tree.xpath("property[@name='sensingstart']/value")[0].text, "%Y-%m-%d %H:%M:%S.%f")
        self.azimuth_line_time = float(tree.xpath("property[@name='azimuthtimeinterval']/value")[0].text)
        self.near_range_slc = float(tree.xpath("property[@name='startingrange']/value")[0].text)
        self.range_pixel_spacing = float(tree.xpath("property[@name='rangepixelsize']/value")[0].text)
        self.azimuth_lines = int(tree.xpath("property[@name='numberoflines']/value")[0].text)
        self.range_samples = int(tree.xpath("property[@name='numberofsamples']/value")[0].text)

def testTxAndSentinel():
    lat = np.arange(-0.1, 0.1, 5e-4) + 40.25
    lon = np.arange(-0.1, 0.1, 5e-4) + 116.3
    h = np.zeros((len(lat), len(lon))) + 50
    lon, lat = np.meshgrid(lon, lat)

    img = Stripmap(
        "/home/sakura/lys/tx/20171206/20171114_slc.xml",
        "/home/sakura/lys/tx/20171206/20171114_slc/20171114.slc",
    )
    y, x = img.geolocate(lon, lat, h)
    stripmap = img.find(y, x)
    stripmap /= np.nanmax(stripmap)

    img = Tops(
        "/home/sakura/lys/20191013/fine_coreg/IW3.xml",
        "/home/sakura/lys/20191013/fine_coreg/IW3/burst_07.slc",
        7
    )
    y, x = img.geolocate(lon, lat, h)
    tops = img.find(y, x)
    tops /= np.nanmax(tops)

    stripmap[270:] = tops[270:]
    plt.imshow(stripmap)
    plt.show()

def testAlosAndSentinel():
    lat = np.arange(-0.15, 0.15, 1e-3) + 40.25
    lon = np.arange(-0.25, 0.25, 1e-3) + 116.4
    # lat = np.arange(-0.1, 0.1, 5e-4) + 40.25
    # lon = np.arange(-0.1, 0.1, 5e-4) + 116.3
    h = np.zeros((len(lat), len(lon))) + 50
    lon, lat = np.meshgrid(lon, lat)

    img = Stripmap(
        "/home/sakura/lys/alos2/20191015/20191015_slc.xml",
        "/home/sakura/lys/alos2/20191015/20191015_slc/20191015.slc",
    )
    y, x = img.geolocate(lon, lat, h)
    stripmap = img.find(y, x)
    stripmap /= np.nanmax(stripmap)

    img = Tops(
        "/home/sakura/lys/20191013/fine_coreg/IW3.xml",
        "/home/sakura/lys/20191013/fine_coreg/IW3/burst_07.slc",
        7
    )
    y, x = img.geolocate(lon, lat, h)
    tops = img.find(y, x)
    tops /= np.nanmax(tops)
    stripmap[160:] = tops[160:]

    plt.imshow(stripmap)
    plt.show()

if __name__ == "__main__":
    # testTxAndSentinel()
    testAlosAndSentinel()