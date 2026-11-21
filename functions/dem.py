from functions.trans import Const, rad2degree, degree2rad
from tifffile import TiffFile
from enum import Enum
import numpy as np
import os


Coord = Enum("Coordinate", ("WGS84", "UTM", "ECS", "Others"))


class Data(object):
    def __init__(self, top: float, left: float, dy: float, dx: float, data: np.ndarray, coord: Coord) -> None:
        """
        注意，数据一定要是第一行对应顶端，最后一行对应底端类型的

        """
        self.top = top
        self.left = left
        self.dy = dy
        self.dx = dx
        self.data = data
        self.height = self.data.shape[0]
        self.width = self.data.shape[1]
        self.coord = coord
        self.__opposite()

    def __opposite(self):
        self.bottom = self.top - self.dy * (self.height - 1)
        self.right = self.left + self.dx * (self.width - 1)

    @property
    def extent(self):
        ys = np.linspace(self.top, self.bottom, self.height)
        xs = np.linspace(self.left, self.right, self.width)
        xs, ys = np.meshgrid(xs, ys)
        return xs, ys

    def find(self, x, y, coord: Coord):
        if self.coord == Coord.WGS84 and coord == Coord.UTM:
            x, y = utm2wgs(x, y)
        elif self.coord == Coord.UTM and coord == Coord.WGS84:
            x, y = wgs2utm(x, y)
        elif self.coord == coord:
            pass
        else:
            raise ValueError("坐标系不匹配")
        i = (self.top - y) / self.dy
        j = (x - self.left) / self.dx
        row, col = np.uint32(i), np.uint32(j)
        rr, cr = i - row, j - col

        rr[row == self.height - 1] = 1
        row[row == self.height - 1] = self.height - 2

        cr[col == self.width - 1] = 1
        col[col == self.width - 1] = self.width - 2

        v1 = (self.data[row, col], self.data[row, col + 1],
              self.data[row + 1, col], self.data[row + 1, col + 1])
        v2 = ((1 - rr) * (1 - cr), (1 - rr) * cr, rr * (1 - cr), rr * cr)
        return np.sum(np.multiply(v1, v2), axis=0, dtype=np.float32)


class Tif(Data):
    def __init__(self, path, coord) -> None:
        with TiffFile(path) as img:
            _, _, _, left, top, _ = img.pages[0].tags.get(
                'ModelTiepointTag').value
            dx, dy, _ = img.pages[0].tags.get('ModelPixelScaleTag').value
            data = img.pages[0].asarray()
            top, left = top - dy / 2, left + dx / 2
            super().__init__(top, left, dy, dx, data, coord)


class Egm(Tif):
    def __init__(self, path) -> None:
        super().__init__(path, Coord.WGS84)
        self.data = np.float32(self.data)

    def find(self, dem: Data):
        x, y = dem.extent
        return super().find(x, y, dem.coord)


class WorldDEM(Tif):
    def __init__(self, path, egm: Egm) -> None:
        super().__init__(path, Coord.WGS84)
        self.data[self.data == -2**15+1] = np.nan
        self.data[:] += egm.find(self)

class DEM(Tif):
    def __init__(self, path, egm: Egm) -> None:
        super().__init__(path, Coord.WGS84)
        self.data = np.float32(self.data)
        self.data[:] += egm.find(self)


class SRTM(Data):
    def __init__(self, path, egm: Egm) -> None:
        dtype = np.dtype(np.int16).newbyteorder('big')
        length = int(np.sqrt(os.path.getsize(path) / 2))
        data = np.fromfile(path, dtype=dtype).reshape((length, length))
        data = np.float32(data)
        d = 1 / (length - 1)
        name = os.path.basename(path)
        top = int(name[1:3]) + 1
        left = int(name[4:7])
        super().__init__(top, left, d, d, data, Coord.WGS84)
        self.data[:] += egm.find(self)


def _center_meridian(zone):
    """
    用于经纬度和UTM投影转换的函数
    """
    return degree2rad(-183 + zone * 6)


def _meridian_length(lat):
    """
    也是用于经纬度和UTM投影转换的函数
    Computes the ellipsoidal distance from the equator to a point at a given latitude.
    Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    """
    n = (Const.a - Const.b) / (Const.a + Const.b)
    alpha = ((Const.a + Const.b) / 2) * (1 + (n ** 2 / 4) + (n ** 4 / 64))
    beta = (-3 * n / 2) + (9 * n ** 3 / 16) + (-3 * n ** 5 / 32)
    gamma = (15 * n ** 2 / 16) + (-15 * n ** 4 / 32)
    delta = (-35 * n ** 3 / 48) + (105 * n ** 5 / 256)
    epsilon = (315 * n ** 4 / 512.0)
    result = alpha * (lat + (beta * np.sin(2 * lat)) + (
        gamma * np.sin(4 * lat)) + (delta * np.sin(
            6 * lat)) + (epsilon * np.sin(8 * lat)))
    return result


def _footpoint_latitude(y):
    """
    还是用于经纬度和UTM投影转换的函数
    """
    alpha_ = ((Const.a + Const.b) / 2) * \
        (1 + (Const.n ** 2 / 4) + (Const.n ** 4 / 64))
    y_ = y / alpha_
    beta_ = (3 * Const.n / 2) + (-27 * Const.n **
                                 3 / 32) + (269 * Const.n ** 5 / 512)
    gamma_ = (21 * Const.n ** 2 / 16) + (-55 * Const.n ** 4 / 32)
    delta_ = (151 * Const.n ** 3 / 96) + (-417 * Const.n ** 5 / 128)
    epsilon_ = (1097 * Const.n ** 4 / 512)
    res = y_ + (beta_ * np.sin(2 * y_)) + (gamma_ * np.sin(4 * y_)) + (delta_ * np.sin(6 * y_)) + (
        epsilon_ * np.sin(8 * y_))
    return res


def wgs2utm(lon, lat):
    """
    WGS84坐标转UTM坐标，注意一定要是角度制表示经纬度

    :param lon: 经度，角度制
    :param lat: 纬度，角度制
    :return: UTM坐标(x, y) 米, ndarray
    """
    zone = 50  # 转成UTM zone 50N
    lon = degree2rad(lon)
    lat = degree2rad(lat)
    center = _center_meridian(zone)
    nu2 = Const.ep2 * np.cos(lat) ** 2
    N = Const.a ** 2 / (Const.b * np.sqrt(1 + nu2))
    t = np.tan(lat)
    t2 = t * t
    l = lon - center
    l3coef = 1 - t2 + nu2
    l4coef = 5 - t2 + 9 * nu2 + 4 * (nu2 * nu2)
    l5coef = 5 - 18 * t2 + (t2 * t2) + 14 * nu2 - 58 * t2 * nu2
    l6coef = 61 - 58 * t2 + (t2 * t2) + 270 * nu2 - 330 * t2 * nu2
    l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2)
    l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2)
    x = N * np.cos(lat) * l + (N / 6.0 * np.cos(lat) ** 3 * l3coef * l ** 3) + (
        N / 120.0 * np.cos(lat) ** 5 * l5coef * l ** 5) + (N / 5040.0 * np.cos(lat) ** 7 * l7coef * l ** 7)

    y = _meridian_length(lat) + (t / 2.0 * N * np.cos(lat) ** 2 * l ** 2) + (
        t / 24 * N * np.cos(lat) ** 4 * l4coef * l ** 4) + (
        t / 720 * N * np.cos(lat) ** 6 * l6coef * l ** 6) + (
        t / 40320 * N * np.cos(lat) ** 8 * l8coef * l ** 8)

    x = x * Const.UTMScaleFactor + 500000
    y = y * Const.UTMScaleFactor
    if y.ndim == 0:
        if y < 0:
            y += 10000000
    else:
        y[y < 0] += 10000000
    return np.array((x, y))


def utm2wgs(x, y):
    """
    UTM坐标转WGS84坐标，注意最终结果以角度制表示经纬度

    :param x: UTM坐标中的x值，米
    :param y: UTM坐标中的y值，米
    :return: 经纬度坐标(lon, lat)，角度制, ndarray
    """
    north, zone = True, 50
    x = (x - 500000.0) / Const.UTMScaleFactor
    if not north:
        y = y - 10000000
    y = y / Const.UTMScaleFactor
    center = _center_meridian(zone)

    lat_f = _footpoint_latitude(y)
    cf = np.cos(lat_f)
    nuf2 = Const.ep2 * cf ** 2
    Nf = Const.a ** 2 / (Const.b * np.sqrt(1 + nuf2))
    Nf2 = Nf * Nf
    tf = np.tan(lat_f)
    tf2 = tf * tf
    tf4 = tf2 * tf2
    x1frac = 1 / (Nf * cf)
    x2frac = tf / (2 * Nf2)
    x3frac = x1frac / (6 * Nf2)
    x4frac = x2frac / (12 * Nf2)
    x5frac = x3frac / (20 * Nf2)
    x6frac = x4frac / (30 * Nf2)
    x7frac = x5frac / (42 * Nf2)
    x8frac = x6frac / (56 * Nf2)
    x2poly = -1 - nuf2
    x3poly = -1 - 2 * tf2 - nuf2
    x4poly = 5 + 3 * tf2 + 6 * nuf2 - 6 * tf2 * nuf2 - \
        3 * (nuf2 * nuf2) - 9 * tf2 * (nuf2 * nuf2)
    x5poly = 5 + 28 * tf2 + 24 * tf4 + 6 * nuf2 + 8 * tf2 * nuf2
    x6poly = -61 - 90 * tf2 - 45 * tf4 - 107 * nuf2 + 162.0 * tf2 * nuf2
    x7poly = -61 - 662 * tf2 - 1320 * tf4 - 720 * (tf4 * tf2)
    x8poly = 1385 + 3633 * tf2 + 4095 * tf4 + 1575 * (tf4 * tf2)
    lat = lat_f + x2frac * x2poly * (
        x * x) + x4frac * x4poly * x ** 4 + x6frac * x6poly * x ** 6 + x8frac * x8poly * x ** 8
    lon = center + x1frac * x + x3frac * x3poly * x ** 3 + \
        x5frac * x5poly * x ** 5 + x7frac * x7poly * x ** 7
    lon = rad2degree(lon)
    lat = rad2degree(lat)
    return np.array((lon, lat))
