from functions.trans import Const, degree2rad, rad2degree, utc2sod, beijing2utc
from tifffile import TiffFile
from lxml import etree
import numpy as np
import re


class Img(object):
    def __init__(self, par, slc_func) -> None:
        self.par = par
        self.__slc = None
        self.__slc_func = slc_func

    def state(self, t, pos=True):
        first_time = self.par['time_of_first_state_vector']
        interval = self.par['state_vector_interval']
        if pos:
            states = np.array(self.par['state_vector_position'])
        else:
            states = np.array(self.par['state_vector_velocity'])
        n_states = 8
        initial = np.int16(np.around((t - first_time) /
                                     interval - (n_states - 1) / 2))
        up, down = len(states) - n_states, 0
        initial[initial < down] = down
        initial[initial > up] = up
        initial = np.expand_dims(initial, axis=0)
        index = np.repeat(initial, n_states, axis=0) + \
            np.arange(n_states).reshape((-1, 1))
        xa = first_time + index * interval
        ya = states[index]
        x = polint(xa, ya[:, :, 0], t)
        y = polint(xa, ya[:, :, 1], t)
        z = polint(xa, ya[:, :, 2], t)
        return np.array((x, y, z))

    def indirect_with_offset(self, lon, lat, h, times=5):
        R_T = wgs2ecs(lon, lat, h).reshape((3, -1))
        t = self.par['start_time'] * np.ones(R_T.shape[1])
        for _ in range(times):
            V_S = self.state(t, False)
            V_S_norm = np.linalg.norm(V_S, axis=0)
            R_ST = R_T - self.state(t, True)
            delta_r = np.sum(np.multiply(R_ST, V_S), axis=0) / V_S_norm
            t += delta_r / V_S_norm

        start_time = self.par['start_time']
        dt = self.par['azimuth_line_time']
        row = (t - start_time) / dt

        pos = self.state(t, True)
        offset = yield pos.reshape((3,) + np.shape(lon))

        target = wgs2ecs(lon, lat, h).reshape((3, -1))
        distance = np.linalg.norm(target-pos, axis=0) + offset.ravel()
        near_range = self.par['near_range_slc']
        dR = self.par['range_pixel_spacing']
        col = (distance - near_range) / dR

        row, col = row.reshape(np.shape(lon)), col.reshape(np.shape(lon))
        if np.max(row) > self.par['azimuth_lines'] - 1:
            row = -2
        elif np.min(row) < 0:
            row = -1
        if np.max(col) > self.par['range_samples'] - 1:
            col = -2
        elif np.min(col) < 0:
            col = -1
        yield row, col
        return

    def indirect(self, lon, lat, h):
        f = self.indirect_with_offset(lon, lat, h)
        pos = f.send(None)
        row, col = f.send(np.zeros(np.shape(lon)))
        return row, col

    def find(self, y, x, threshold=5):
        if np.any(y < 0) or np.any(x < 0):
            raise ValueError("超出影像边界")
        if self.__slc is None:
            self.__slc = self.__slc_func()
        height, width = self.__slc.shape

        row, col = np.uint32(y), np.uint32(x)
        rr, cr = y - row, x - col
        rr[row == height - 1] = 1
        row[row == height - 1] = height - 2
        cr[col == width - 1] = 1
        col[col == width - 1] = width - 2
        v1 = (self.__slc[row, col], self.__slc[row, col + 1],
              self.__slc[row + 1, col], self.__slc[row + 1, col + 1])
        v2 = ((1 - rr) * (1 - cr), (1 - rr) * cr, rr * (1 - cr), rr * cr)
        intensity = np.sum(np.multiply(np.abs(v1), v2),
                           axis=0, dtype=np.float32)

        max_v = np.nanpercentile(intensity, 100-threshold)
        min_v = np.nanpercentile(intensity, threshold)
        intensity[intensity > max_v] = max_v
        intensity[intensity < min_v] = min_v
        return intensity
    
    def set(self, row, col, color=10000):
        if self.__slc is None:
            self.__slc = self.__slc_func()
        self.__slc[int(row), int(col)] = color

class TX(Img):
    def __init__(self, xml_path, cos_path) -> None:
        par_str = tx_xml2par(xml_path)
        par = read_par(par_str)
        super().__init__(par, lambda: tx_cos2slc(cos_path))
        self.f = 9.65e9


class TH(Img):
    def __init__(self, xml_path, cos_path) -> None:
        par_str = th_xml2par(xml_path)
        par = read_par(par_str)
        super().__init__(par, lambda: th_cos2slc(cos_path))
        self.f = 9.65e9

    def indirect_with_offset(self, lon, lat, h, ref=None, times=10):
        if ref is None:
            ref = self
        c = 2.99792458e8
        R_T = wgs2ecs(lon, lat, h).reshape((3, -1))
        t = self.par['start_time'] * np.ones(R_T.shape[1])
        def func_Rrt(t):return R_T - ref.state(t, True)
        def func_Rst(t):return R_T- self.state(t,True)
        def func_norm(m): return np.linalg.norm(m,axis=0)
        def func_dot(m1, m2): return np.sum(np.multiply(m1, m2), axis=0)
        for _ in range(times):
            Vr = ref.state(t - func_norm(func_Rrt(t)) / c, False)
            Vs = self.state(t + func_norm(func_Rst(t)) / c,False)

            t += (func_dot(func_Rrt(t - func_norm(func_Rrt(t)) / c), Vr)/func_norm(Vr)**2 + func_dot(func_Rst(t + func_norm(func_Rst(t)) / c), Vs)/func_norm(Vs)**2) / 2

        start_time = self.par['start_time']
        dt = self.par['azimuth_line_time']
        row = (t - start_time + func_norm(func_Rst(t)) / c) / dt

        offset = yield self.state(t, True).reshape((3,) + np.shape(lon))
        dR = self.par['range_pixel_spacing']
        
        tr = t - func_norm(func_Rrt(t)) / c
        ts = t + func_norm(func_Rst(t)) / c
        col = ((func_norm(func_Rrt(tr))+func_norm(func_Rst(ts)))/2 - self.par['near_range_slc'] + offset.ravel()) / dR


        row, col = row.reshape(np.shape(lon)), col.reshape(np.shape(lon))
        if np.max(row) > self.par['azimuth_lines'] - 1:
            row = -2
        elif np.min(row) < 0:
            row = -1
        if np.max(col) > self.par['range_samples'] - 1:
            col = -2
        elif np.min(col) < 0:
            col = -1

        yield row, col
        return

class S1(Img):
    xml_path = None
    tif = None

    def __init__(self, xml_path, eof_path, tif_path, burst) -> None:
        par = s1_xml2par(xml_path, eof_path, burst)

        def slc_func():
            if xml_path != S1.xml_path:
                S1.xml_path = xml_path
                S1.tif = s1_tif2slc(tif_path)
            return S1.tif[par['line_offset']:par['line_offset']+par['azimuth_lines']]
        super().__init__(par, slc_func)
        self.f = 5.41e9


class Alos(Img):
    def __init__(self, par_path, slc_path) -> None:
        with open(par_path, 'r') as f:
            par = read_par(f.read())
            super().__init__(par, lambda: read_slc(slc_path, par))
            self.f = 1.26e9


def wgs2ecs(lon, lat, height):
    """
    经纬度转地心直角坐标

    :param lon: 经度，角度制
    :param lat: 纬度，角度制
    :param height: 高程，米
    :return: [x, y, z] (ndarray, 单位 米)
    """
    lon = degree2rad(lon)
    lat = degree2rad(lat)
    N = Const.a / np.sqrt(1 - Const.e2 * np.sin(lat) ** 2)
    x = (N + height) * np.cos(lat) * np.cos(lon)
    y = (N + height) * np.cos(lat) * np.sin(lon)
    z = (N * (1 - Const.e2) + height) * np.sin(lat)
    return np.array((x, y, z))


def ecs2wgs(x, y, z):
    """
    地心直角坐标转经纬度

    :param x: x坐标值，米
    :param y: y坐标值，米
    :param z: z坐标值，米
    :return: [lon, lat, height] （经度，纬度，高程）（角度制，米）
    """
    p = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan(z * Const.a / (p * Const.b))
    lon = np.arctan2(y, x)
    lat = np.arctan((z + Const.e2 / (1 - Const.e2) * Const.b * np.sin(theta) ** 3
                     ) / (p - Const.e2 * Const.a * np.cos(theta) ** 3))
    N = Const.a / np.sqrt(1 - Const.e2 * np.sin(lat) ** 2)
    height = p / np.cos(lat) - N
    lon, lat = rad2degree(lon), rad2degree(lat)
    return np.array((lon, lat, height))


def read_par(par: str):
    """
    读取par的内容，从头文件中提取出键对应的值

    :param par: par文件的内容，或者文件路径
    :return: 包含内容的字典
    """
    keys = ['start_time', 'azimuth_line_time', 'near_range_slc', 'range_pixel_spacing',
            'state_vector_position', 'state_vector_velocity', 'time_of_first_state_vector',
            'state_vector_interval']
    keys += ['range_samples', 'azimuth_lines', 'image_format']

    def str2num(s):
        """
        若字符串实际为数字，则转换为数字；否则不变。
        能力有限，函数通用性不强，针对gamma头文件格式特化。

        :param s: 字符串
        :return: 数字或字符串
        """

        if '.' in s:
            return float(s)
        elif re.fullmatch(r'[0-9]+', s):
            return int(s)
        else:
            return s

    res = {}
    for i in par.split("\n"):  # 遍历所有头文件的参数
        for k in keys:  # 遍历所有所需键
            if i[:len(k)] == k:  # 输入的键只要能和参数全名的前一部分相同就行了，可能会省略编号
                ans = i.split()[1:]  # 按空格切分
                if len(ans) == 1 or len(ans) == 2:
                    # 该情况表明ans是个单词/数值（1），或者带单位（2），但总之不是多个数值。不存在多个单词的情况
                    ans = str2num(ans[0])
                else:  # 否则是多个数值，该情况下必定带有单位，因此取前半部分就行了
                    ans = [str2num(ans[i])
                           for i in range(int(len(ans) / 2))]
                if k not in res.keys():  # 如果之前没有添加过该key。一般情况下都是这种情况。
                    # 例外是有多组值的key（比如state_vector_position），第2次及以后的提取。
                    res[k] = ans
                else:  # 对于多组值第2次及以后的情况
                    if not isinstance(res[k][0], list):  # 第2次加入，键值不是列表而是第1次提取的值
                        res[k] = [res[k], ans]  # 将第1次的值放入新建的列表
                    else:  # 第3次及以后加入
                        res[k].append(ans)  # 直接加入列表
    if len(res) != len(keys):  # 如果没能把所有键都提取出来，也不能哑播放，要报错
        k = []
        for i in keys:
            if i not in res.keys():
                k.append(i)
        raise ValueError("没有这些量：%s" % ",".join(k))
    return res


def read_slc(slc_path, par_dict):
    """
    给定slc路径，读取slc数据为numpy数组。千万注意Gamma是Linux软件，二进制文件为大端序。

    :params slc_path: slc文件路径
    :params par_dict: par的解析结果
    :return: ndarray
    """
    cols, rows, dtype = par_dict['range_samples'], par_dict['azimuth_lines'], par_dict['image_format']
    if dtype == 'SCOMPLEX':
        dtype = np.dtype(np.int16).newbyteorder('big')
    elif dtype == 'FCOMPLEX':
        dtype = np.dtype(np.float32).newbyteorder('big')
    else:
        raise Exception('未知类型的image_format')
    res = np.float32(np.fromfile(
        slc_path, dtype=dtype).reshape((rows, cols * 2)))
    return res[:, ::2] + res[:, 1::2]*1j


def polint(xa, ya, x):
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


def tx_xml2par(xml_path):
    """
    将xml元数据文件转换gamma采用的元数据格式
    """
    par = {key: '' for key in get_keys(formatted=False)}
    tree = etree.parse(xml_path)
    root = tree.getroot()
    productInfo = root.find('productInfo')
    sceneInfo = productInfo.find('sceneInfo')
    par['title'] = sceneInfo.find('sceneID').text
    par['sensor'] = productInfo.find('missionInfo').find('mission').text
    start_UTC = sceneInfo.find('start').find('timeUTC').text
    date = re.split('[-T:Z]', start_UTC)[:6]
    par['date'] = '%d %d %d %d %d %6.4f' % (int(date[0]), int(date[1]), int(
        date[2]), int(date[3]), int(date[4]), float(date[5]))
    start_time = utc2sod(start_UTC)
    par['start_time'] = '%15.6f   s' % start_time
    end_time = utc2sod(sceneInfo.find('stop').find('timeUTC').text)
    if end_time < start_time:
        end_time += 24 * 3600
    center_time = (start_time + end_time) / 2
    par['center_time'] = '%15.6f   s' % center_time
    par['end_time'] = '%15.6f   s' % end_time
    complexImageInfo = root.find(
        'productSpecific').find('complexImageInfo')
    prf = float(complexImageInfo.find('commonPRF').text)
    par['azimuth_line_time'] = '%15.7e   s' % (1 / prf)
    par['line_header_size'] = '%12d' % 0  # 源代码如此
    imageDataInfo = productInfo.find('imageDataInfo')
    imageRaster = imageDataInfo.find('imageRaster')
    par['range_samples'] = '%12d' % int(
        imageRaster.find('numberOfColumns').text)
    par['azimuth_lines'] = '%12d' % int(
        imageRaster.find('numberOfRows').text)
    par['range_looks'] = '%12d' % float(
        imageRaster.find('rangeLooks').text)
    par['azimuth_looks'] = '%12d' % float(
        imageRaster.find('azimuthLooks').text)
    imageDataType = imageDataInfo.find('imageDataType').text
    imageDataFormat = imageDataInfo.find('imageDataFormat').text
    if imageDataType == 'COMPLEX' and imageDataFormat == 'COSAR':
        par['image_format'] = 'SCOMPLEX'
    else:
        raise Exception('ERROR: unknown image data type and format: %s %s' % (
            imageDataType, imageDataFormat))
    projection = productInfo.find(
        'productVariantInfo').find('projection').text
    if projection == 'SLANTRANGE':
        par['image_geometry'] = 'SLANT_RANGE'
    else:
        raise Exception('ERROR: unknown image projection: %s' % projection)
    par['range_scale_factor'] = '%15.7e' % 1.0  # 源代码如此
    par['azimuth_scale_factor'] = '%15.7e' % 1.0  # 源代码如此
    sceneCenterCoord = sceneInfo.find('sceneCenterCoord')
    center_latitude = float(sceneCenterCoord.find('lat').text)
    center_longitude = float(sceneCenterCoord.find('lon').text)
    par['center_latitude'] = '%14.7f   degrees' % center_latitude
    par['center_longitude'] = '%14.7f   degrees' % center_longitude
    par['heading'] = '%14.7f   degrees' % float(
        sceneInfo.find('headingAngle').text)
    par['range_pixel_spacing'] = '%12.6f   m' % float(
        complexImageInfo.find('projectedSpacingRange').find('slantRange').text)
    par['azimuth_pixel_spacing'] = '%12.6f   m' % float(
        complexImageInfo.find('projectedSpacingAzimuth').text)
    rangeTime = sceneInfo.find('rangeTime')
    c = 2.99792458e8
    near_range_slc = float(rangeTime.find('firstPixel').text) * c / 2
    far_range_slc = float(rangeTime.find('lastPixel').text) * c / 2
    par['near_range_slc'] = '%15.4f  m' % near_range_slc
    par['center_range_slc'] = '%15.4f  m' % (
        (near_range_slc + far_range_slc) / 2)
    par['far_range_slc'] = '%15.4f  m' % far_range_slc
    par['first_slant_range_polynomial'] = '%12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3 ' % (
        0, 0, 0, 0, 0, 0)  # 源代码如此
    par['center_slant_range_polynomial'] = '%12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3 ' % (
        0, 0, 0, 0, 0, 0)  # 源代码如此
    par['last_slant_range_polynomial'] = '%12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3 ' % (
        0, 0, 0, 0, 0, 0)  # 源代码如此
    par['incidence_angle'] = '%12.4f   degrees' % float(
        sceneCenterCoord.find('incidenceAngle').text)
    par['azimuth_deskew'] = 'ON'  # 源代码如此
    par['azimuth_angle'] = '%12.4f   degrees' % (90 if productInfo.find(
        'acquisitionInfo').find('lookDirection').text == 'RIGHT' else -90)
    instrument = root.find('instrument')
    par['radar_frequency'] = '%15.7e  Hz' % float(
        instrument.find('radarParameters').find('centerFrequency').text)
    par['adc_sampling_rate'] = '%15.7e  Hz' % float(
        complexImageInfo.find('commonRSF').text)
    processingParameter = root.find(
        'processing').find('processingParameter')
    par['chirp_bandwidth'] = '%15.7e  Hz' % float(
        processingParameter.find('totalProcessedRangeBandwidth').text)
    par['prf'] = '%12.7f  Hz' % prf
    par['azimuth_proc_bandwidth'] = '%12.5f  Hz' % float(
        processingParameter.find('totalProcessedAzimuthBandwidth').text)
    averageDopplerCentroidInfo = root.find('productQuality').find(
        'processingParameterQuality').find('averageDopplerCentroidInfo')
    azimuthCenter = averageDopplerCentroidInfo.find('azimuthCenter')
    dp0 = float(azimuthCenter.find('midRange').text)
    dp1 = (float(azimuthCenter.find('farRange').text) -
           float(azimuthCenter.find('nearRange').text)) / (far_range_slc - near_range_slc)
    par['doppler_polynomial'] = '%12.5f %12.5e %12.5e %12.5e  Hz     Hz/m     Hz/m^2     Hz/m^3' % (
        dp0, dp1, 0, 0)  # 源代码如此
    par['doppler_poly_dot'] = '%12.5e %12.5e %12.5e %12.5e  Hz/s   Hz/s/m   Hz/s/m^2   Hz/s/m^3' % (
        (float(averageDopplerCentroidInfo.find('azimuthStop').find('midRange').text) - float(
            averageDopplerCentroidInfo.find('azimuthStart').find('midRange').text)) / (
                end_time - start_time), 0, 0, 0)  # 源代码如此
    par['doppler_poly_ddot'] = '%12.5e %12.5e %12.5e %12.5e  Hz/s^2 Hz/s^2/m Hz/s^2/m^2 Hz/s^2/m^3' % (
        0, 0, 0, 0)  # 源代码如此
    par['receiver_gain'] = '%12.4f  dB' % abs(10 * np.log10(float(instrument.find(
        'settings').find('rxGainSetting').find('rxGain').text)))
    par['calibration_gain'] = '%12.4f  dB' % abs(10 * np.log10(float(root.find(
        'calibration').find('calibrationConstant').find('calFactor').text)))
    orbit = root.find('platform').find('orbit')
    orbitHeader = orbit.find('orbitHeader')
    number_of_state_vectors = int(orbitHeader.find('numStateVectors').text)
    time_of_first_state_vector = utc2sod(orbitHeader.find(
        'firstStateTime').find('firstStateTimeUTC').text)
    state_vector_interval = float(
        orbitHeader.find('stateVectorTimeSpacing').text)
    state_vectors = [{i: float(node.find(i).text) for i in [
        'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']} for node in orbit.findall('stateVec')]
    par['sar_to_earth_center'] = '%15.4f   m' % np.linalg.norm(tx_interp(
        center_time, time_of_first_state_vector, state_vector_interval, state_vectors))
    earth_semi_major_axis = 6378137.0  # 源代码如此
    earth_semi_minor_axis = 6356752.3141  # 源代码如此

    TERRA_ALT = 200.0  # /* nominal terrain altitude */
    par['earth_radius_below_sensor'] = '%15.4f   m' % np.linalg.norm(
        wgs2ecs(center_latitude, center_longitude, TERRA_ALT))
    par['earth_semi_major_axis'] = '%15.4f   m' % earth_semi_major_axis
    par['earth_semi_minor_axis'] = '%15.4f   m' % earth_semi_minor_axis
    par['number_of_state_vectors'] = '%15d' % number_of_state_vectors
    # 以本数据为例，理应是12，但gamma软件的结果是11，state_vectors也只有11个，但应该不影响最终结果
    par['time_of_first_state_vector'] = '%15.6f   s' % time_of_first_state_vector
    par['state_vector_interval'] = '%15.6f   s' % state_vector_interval

    keys = get_keys(formatted=True)
    res = ""
    res += 'Gamma Interferometric SAR Processor (ISP) - Image Parameter File\n\n'
    for k, v in zip(keys, par.values()):
        res += ('%s%s\n' % (k, v))
    for i, j in enumerate(state_vectors):
        res += ('state_vector_position_%d: %14.4f  %14.4f  %14.4f   m   m   m\n' % (
            i + 1, j['posX'], j['posY'], j['posZ']))
        res += ('state_vector_velocity_%d: %14.5f  %14.5f  %14.5f   m/s m/s m/s\n' % (
            i + 1, j['velX'], j['velY'], j['velZ']))
    return res


def get_keys(formatted):
    keys = [
        'title:     ',
        # ascii string with title of the scene, e.g., C85_N22_A_SM_stripNear_011_R_2009-12-26T12:36:00.293000Z
        'sensor:    ',
        # sensor name (RADARSAT, SRL-1, SRL-2, ERS-1, ERS-2, JERS-1,...), e.g., TSX-1
        'date:      ',
        # date in form: YYYY MM DD hh mm ss.ttt UTC, e.g., 2009 12 26 12 36 0.2930
        'start_time:          ',
        # time of image start UTC seconds since start of day, e.g., 45360.293000   s
        'center_time:         ',
        # time of image center UTC seconds since start of day, e.g., 45364.292979   s
        'end_time:            ',
        # time of image end UTC seconds since start of day, e.g., 45368.292959   s
        'azimuth_line_time:   ',
        # time per azimuth line (s), e.g., 3.2324374e-04   s
        'line_header_size:       ',
        # header size in bytes for each image line, e.g., 0
        'range_samples:          ',
        # number of range pixels/line in the image, e.g., 11342
        'azimuth_lines:          ',
        # number of range lines in the scene, e.g., 24750
        'range_looks:            ',
        # number of range looks, for SLC = 1, e.g., 1
        'azimuth_looks:          ',
        # number of azimuth looks, for SLC = 1, e.g., 1
        'image_format:               ',
        # image format FCOMPLEX, SCOMPLEX, FLOAT, SHORT, BYTE, e.g., SCOMPLEX
        'image_geometry:             ',
        # image geometry type SLANT_RANGE, GROUND_RANGE, GEOCODED, e.g., SLANT_RANGE
        'range_scale_factor:   ',
        # range pixel spacing scale factor, without resampling rpssf=1.0, e.g., 1.0000000e+00
        'azimuth_scale_factor: ',
        # azimuth pixel spacing scale factor, without resampling azpssf=1.0, e.g., 1.0000000e+00
        'center_latitude:      ',
        # latitude of scene center in decimal degrees, e.g., 28.3242204   degrees
        'center_longitude:     ',
        # longitude of scene center in decimal degrees, e.g., 80.8263537   degrees
        'heading:              ',
        # sub-satellite track heading at scene center (decimal-degrees), e.g., 349.9712523   degrees
        'range_pixel_spacing:    ',
        # slant range pixel spacing  (meters), e.g., 0.909403   m
        'azimuth_pixel_spacing:  ',
        # azimuth along track pixel spacing (meters), e.g., 2.288608   m
        'near_range_slc:       ',
        # near slant of image (meters), e.g., 635702.2645  m
        'center_range_slc:     ',
        # center slant of image (meters), e.g., 640859.0369  m
        'far_range_slc:        ',
        # far slant of image (meters), e.g., 646015.8093  m
        'first_slant_range_polynomial:   ',
        # first slant range polynomial coefficients. 可能在slc影像中是没用的，在地距影像中使用
        # 令 first_slant_range_polynomial = f[]
        # f[0] contains the reference orbit time for the polynomial
        # slant range = f[1] + f[2]*(GR-r0) + f[3]*(GR-r0)**2 + f[4]*(GR-r0)**3 + f[5]*(GR-r0)**4
        # (r0 is the ground range of the first pixel)
        # e.g., 0.00000      0.00000  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  s m 1 m^-1 m^-2 m^-3
        'center_slant_range_polynomial:  ',
        # center slant range polynomial coefficients
        # e.g., 0.00000      0.00000  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  s m 1 m^-1 m^-2 m^-3
        'last_slant_range_polynomial:    ',
        # last slant range polynomial coefficients
        # e.g., 0.00000      0.00000  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  s m 1 m^-1 m^-2 m^-3
        'incidence_angle:        ',
        # incidence angle at the center of the scene (deg.), e.g., 38.7401   degrees
        'azimuth_deskew:          ',
        # azimuth deskew (ON, OFF), e.g., ON (零多普勒化)
        'azimuth_angle:          ',
        # nominal azimuth antenna angle (decimal degrees CW about N, right looking SAR 90.0, left looking: -90.0)
        # e.g., 90.0000   degrees
        'radar_frequency:      ',
        # radar carrier center frequency (Hz), e.g., 9.6499993e+09  Hz
        'adc_sampling_rate:    ',
        # sample rate  of radar analog to digital converter (Hz), e.g., 1.6482918e+08  Hz
        'chirp_bandwidth:      ',
        # radar range chirp bandwidth (Hz), e.g., 1.5000000e+08  Hz
        'prf:                     ',
        # radar pulse repetition frequency (Hz), e.g., 3093.6407742  Hz
        'azimuth_proc_bandwidth:  ',
        # 3 dB azimuth processing bandwidth (Hz), e.g., 1380.00000  Hz
        'doppler_polynomial:     ',
        # doppler centroid polynomial coefficients. 令 doppler_polynomial = dp[]
        # dp[0] + dp[1]*(r-r1) + dp[2]*(r-r1)**2 + dp[3]*(r-r1)**3
        # (r is the slant range and r1 is the slant range at the center of the image
        # e.g., -133.02999 -9.87748e-06  0.00000e+00  0.00000e+00  Hz     Hz/m     Hz/m^2     Hz/m^3
        'doppler_poly_dot:       ',
        # derivative w.r.t. along-track time of each of the terms in doppler_polynomial
        # e.g., -7.72281e-01  0.00000e+00  0.00000e+00  0.00000e+00  Hz/s   Hz/s/m   Hz/s/m^2   Hz/s/m^3
        'doppler_poly_ddot:      ',
        # second derivative w.r.t. along-track time of each of the terms in doppler_polynomial
        # e.g., 0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  Hz/s^2 Hz/s^2/m Hz/s^2/m^2 Hz/s^2/m^3
        'receiver_gain:           ',
        # receiver gain (dB), e.g., 23.7770  dB
        'calibration_gain:        ',
        # calibration gain (dB), e.g., 50.0623  dB
        'sar_to_earth_center:          ',
        # distance of SAR sensor from earth center at center scene, e.g., 6885140.9799   m
        'earth_radius_below_sensor:    ',
        # center scene geocentric radius (m), e.g., 6373554.2808   m
        'earth_semi_major_axis:        ',
        # earth ellipsoid semi-major axises (m), e.g., 6378137.0000   m
        'earth_semi_minor_axis:        ',
        # earth ellipsoid semi-minor axises (m), e.g., 6356752.3141   m
        'number_of_state_vectors:      ',
        # number of state vectors, e.g., 11
        'time_of_first_state_vector:   ',
        # UTC time (sec) since start of day for first state vector, e.g., 45315.000000   s
        'state_vector_interval:        ',
        # time interval between state vectors (s), e.g., 10.000000   s
        # 最后写入文件时还需要加入 state vector, e.g.,
        # state_vector_position_1:   1321998.5639    6126956.7235    2850514.1229   m   m   m
        # state_vector_velocity_1:      845.13700     -3383.02300      6852.76700   m/s m/s m/s
    ]
    return keys if formatted else [key.split(':')[0] for key in keys]


def tx_interp(center_time, time_of_first_state_vector, state_vector_interval, state_vectors):
    """
    卫星状态插值，用来计算'sar_to_earth_center'，center_time时的卫星到地心距离
    这里看不懂没关系，因为参数对xml文件特化了。定位算法中有通用的轨道插值函数
    :param center_time: 意义同par文件中的同名参数
    :param time_of_first_state_vector: 意义同par文件中的同名参数
    :param state_vector_interval: 意义同par文件中的同名参数
    :param state_vectors: 意义同par文件中的同名参数
    :return: center_time时刻的卫星位置(x,y,z)
    """
    number_of_state_vectors = len(state_vectors)
    center_first = center_time - time_of_first_state_vector
    center_last = center_first - \
        (number_of_state_vectors - 1) * state_vector_interval
    MAX_INTERP_TIME = 90.0  # range of orbit propagation +/- 90.0 sec
    if center_first < -MAX_INTERP_TIME:
        raise Exception(
            'ERROR: center_time less than range of interpolation (s): %12.5f' % center_time)
    elif center_last > MAX_INTERP_TIME:
        raise Exception(
            'ERROR: center_time greater than range of interpolation (s): %12.5f' % center_time)
    ppts = np.min([8, number_of_state_vectors])
    # number of points used for interpolation
    initial = np.around(
        center_first / state_vector_interval - (ppts - 1) / 2).astype(np.int32)
    max_initial = number_of_state_vectors - ppts
    if initial < 0:
        initial = 0
    elif initial > max_initial:
        initial = max_initial
    interval = range(initial, initial + ppts)
    pos = {i: np.array([state_vectors[j][i] for j in interval]).reshape((-1, 1))
           for i in ('posX', 'posY', 'posZ')}
    t = np.array([time_of_first_state_vector + i *
                  state_vector_interval for i in interval]).reshape((-1, 1))
    c = np.array([1]) * center_time
    return (polint(t, pos['posX'], center_time), polint(
        t, pos['posY'], c), polint(t, pos['posZ'], c))


def tx_cos2slc(cos_path):
    endian = 'big'  # TerraSAR的数据是大端序存储的
    with open(cos_path, 'rb') as cos:
        s = 4
        Bytes_In_Burst = int.from_bytes(
            cos.read(s), byteorder=endian)  # BIB, 一个burst的总字节数
        Range_Sample_Relative_Index = int.from_bytes(
            cos.read(s), byteorder=endian)
        # RSRI, 距离向采样点索引值，给出了在一个虚拟栅格上每一行距离向上第一个采样点相对于图像中心参考点零多普勒时刻位置的索引值
        Range_Samples = int.from_bytes(cos.read(s), byteorder=endian)
        # RS, 距离向采样点数，该值在所有burst中保持一致
        Azimuth_Samples = int.from_bytes(
            cos.read(s), byteorder=endian)  # AS, 方位向采样点数
        Burst_Index = int.from_bytes(
            cos.read(s), byteorder=endian)  # BI, burst个数的索引值
        Rangeline_Total_Number_of_Bytes = int.from_bytes(
            cos.read(s), byteorder=endian)  # RTNB, 一行距离向的总字节数，包含注释信息
        Total_Number_of_Lines = int.from_bytes(cos.read(s), byteorder=endian)
        # TNL, 方位向上的长度，即行数，包含注释信息。该值只在文件第一行给出，方便文件读取和其他burst替换。因此整个cos文件大小可以通过RTNB×TNL得出
        # print('Bytes_In_Burst:', Bytes_In_Burst)
        # print('Range_Sample_Relative_Index:', Range_Sample_Relative_Index)
        # print('Range_Samples:', Range_Samples)
        # print('Azimuth_Samples:', Azimuth_Samples)
        # print('Burst_Index:', Burst_Index)
        # print('Rangeline_Total_Number_of_Bytes:',
        #       Rangeline_Total_Number_of_Bytes)
        # print('Total_Number_of_Lines:', Total_Number_of_Lines)
        # 这里Burst_Index = 1，条带模式的图像数据全都是一个burst，我也不知道不为1时是什么样子的
        annotation_lines = Total_Number_of_Lines - Azimuth_Samples
        annotation_columns = Rangeline_Total_Number_of_Bytes // s - Range_Samples
        cos.seek(annotation_lines * Rangeline_Total_Number_of_Bytes)
        dtype = np.dtype(np.int16).newbyteorder(endian)
        res = np.float32(np.frombuffer(cos.read(Rangeline_Total_Number_of_Bytes * Azimuth_Samples), dtype=dtype).reshape(
            Azimuth_Samples, (Range_Samples + annotation_columns)*2)[:, 4:])
        return res[:, ::2]+res[:, 1::2]*1j


def th_xml2par(xml_path):
    """
    将天绘二号卫星影像的头文件.xml转换为gamma通用的par头文件

    :param xml_path: .xml文件的读取路径
    :return: 无，结果保存至par头文件字符串
    """
    par = {key: '' for key in get_keys(False)}  # 去掉多余的冒号及之后的空格
    tree = etree.parse(xml_path)  # 将xml头文件转换为etree，而不是字符串，以便索引
    product = tree.getroot()  # 根节点的标签就是product
    par['title'] = product.find('sceneID').text
    par['sensor'] = product.find('satellite').text
    imageInfo = product.find('imageInfo')
    imageTime = imageInfo.find('imageTime')
    start_UTC = beijing2utc(imageTime.find('start').text)
    date = re.split('[-T:]', start_UTC)
    par['date'] = '%d %d %d %d %d %6.4f' % (int(date[0]), int(date[1]), int(
        date[2]), int(date[3]), int(date[4]), float(date[5]))
    start_time = utc2sod(start_UTC)
    par['start_time'] = '%15.6f   s' % start_time
    end_time = utc2sod(beijing2utc(imageTime.find('end').text))
    if end_time < start_time:
        end_time += 24 * 3600
    center_time = (start_time + end_time) / 2
    par['center_time'] = '%15.6f   s' % center_time
    par['end_time'] = '%15.6f   s' % end_time
    sensor = product.find('sensor')
    waveParams = sensor.find('waveParams')
    prf = float(waveParams.find('prf').text)
    par['azimuth_line_time'] = '%15.7e   s' % (1 / prf)
    par['line_header_size'] = '%12d' % 0
    par['range_samples'] = '%12d' % int(imageInfo.find('width').text)
    par['azimuth_lines'] = '%12d' % int(imageInfo.find('height').text)
    processInfo = product.find('processInfo')
    par['range_looks'] = '%12d' % int(
        processInfo.find('MultilookRange').text)
    par['azimuth_looks'] = '%12d' % int(
        processInfo.find('MultilookAzimuth').text)
    par['image_format'] = 'SCOMPLEX'
    par['image_geometry'] = 'SLANT_RANGE'
    par['range_scale_factor'] = '%15.7e' % 1.0
    par['azimuth_scale_factor'] = '%15.7e' % 1.0
    center = imageInfo.find('center')
    center_latitude = float(center.find('latitude').text)
    center_longitude = float(center.find('longitude').text)
    par['center_latitude'] = '%14.7f   degrees' % center_latitude
    par['center_longitude'] = '%14.7f   degrees' % center_longitude
    platform = product.find('platform')
    par['heading'] = '%14.7f   degrees' % float(
        platform.find('YawAngle').text)
    par['range_pixel_spacing'] = '%12.6f   m' % float(
        imageInfo.find('widthSpace').text)
    par['azimuth_pixel_spacing'] = '%12.6f   m' % float(
        imageInfo.find('heightSpace').text)
    near_range_slc = float(imageInfo.find('nearRange').text)
    center_range_slc = float(imageInfo.find('refRange').text)
    delta_range_slc = center_range_slc - near_range_slc
    far_range_slc = center_range_slc + delta_range_slc
    par['near_range_slc'] = '%15.4f  m' % near_range_slc
    par['center_range_slc'] = '%15.4f  m' % center_range_slc
    par['far_range_slc'] = '%15.4f  m' % far_range_slc
    par['first_slant_range_polynomial'] = '%12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3 ' % (
        0, 0, 0, 0, 0, 0)
    par['center_slant_range_polynomial'] = '%12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3 ' % (
        0, 0, 0, 0, 0, 0)
    par['last_slant_range_polynomial'] = '%12.5f %12.5f %12.5e %12.5e %12.5e %12.5e  s m 1 m^-1 m^-2 m^-3 ' % (
        0, 0, 0, 0, 0, 0)
    par['incidence_angle'] = '%12.4f   degrees' % float(
        waveParams.find('centerAngle').text)
    par['azimuth_deskew'] = 'ON'
    par['azimuth_angle'] = '%12.4f   degrees' % (
        90 if sensor.find('lookDirection').text == 'R' else -90)
    par['radar_frequency'] = '%15.7e  Hz' % (float(
        product.find('carrierFrequency').text) * 1e9)
    par['adc_sampling_rate'] = '%15.7e  Hz' % (float(
        waveParams.find('sampleRate').text) * 1e6)
    par['chirp_bandwidth'] = '%15.7e  Hz' % (float(
        waveParams.find('bandWidth').text) * 1e6)
    par['prf'] = '%12.7f  Hz' % prf
    par['azimuth_proc_bandwidth'] = '%12.5f  Hz' % float(
        waveParams.find('azimuthBandwidth').text)
    AzFdc0 = float(processInfo.find('AzFdc0').text)
    # 在TerraSAR的数据中，下几行的那个d0参数就是center_range_slc处的多普勒频率
    # 但是，天绘数据给出了一个独立的参考快时间，即AzFdc0，这不是center_range_slc对应的快时间
    # 因此，d0需要重新计算，把参考斜距调整到center_range_slc
    DopplerCentroidCoefficients = processInfo.find(
        'DopplerCentroidCoefficients')
    c = 2.99792458e8
    d0 = float(DopplerCentroidCoefficients.find('d0').text) + float(
        DopplerCentroidCoefficients.find('d1').text) * (center_range_slc * 2 / c - AzFdc0)
    d1 = float(DopplerCentroidCoefficients.find('d1').text) * 2 / c
    par['doppler_polynomial'] = '%12.5f %12.5e %12.5e %12.5e  Hz     Hz/m     Hz/m^2     Hz/m^3' % (
        d0, d1, 0, 0)
    par['doppler_poly_dot'] = '%12.5e %12.5e %12.5e %12.5e  Hz/s   Hz/s/m   Hz/s/m^2   Hz/s/m^3' % (
        0, 0, 0, 0)
    # doppler_poly_dot是对最终产品而言的多普勒频率变化率。举个例子，在TerraSAR中，是利用Azimuth行号最小和最大处
    # 两行的center_range处的多普勒频率相减，然后除以成像时间，即6.2Hz/8s=0.77Hz/s的样子，系数很小。
    # 但是天绘卫星给出来的是卫星运动造成的多普勒调频率，数量级是10^3Hz/s。因此这个数据是缺失的。
    par['doppler_poly_ddot'] = '%12.5e %12.5e %12.5e %12.5e  Hz/s^2 Hz/s^2/m Hz/s^2/m^2 Hz/s^2/m^3' % (
        0, 0, 0, 0)
    par['receiver_gain'] = '%12.4f  dB' % float(
        waveParams.find('receiverGain').text)
    par['calibration_gain'] = '%12.4f  dB' % 1
    par['sar_to_earth_center'] = '%15.4f   m' % float(
        platform.find('Rs').text)
    earth_semi_major_axis = 6378137.0
    earth_semi_minor_axis = 6356752.3141

    TERRA_ALT = 200.0  # nominal terrain altitude
    par['earth_radius_below_sensor'] = '%15.4f   m' % np.linalg.norm(
        wgs2ecs(center_latitude, center_longitude, TERRA_ALT))
    par['earth_semi_major_axis'] = '%15.4f   m' % earth_semi_major_axis
    par['earth_semi_minor_axis'] = '%15.4f   m' % earth_semi_minor_axis
    GPSParam = product.find('GPS').findall('GPSParam')
    par['number_of_state_vectors'] = '%15d' % len(GPSParam)
    time_of_first_state_vector = utc2sod(
        beijing2utc(GPSParam[0].find('TimeStamp').text))
    par['time_of_first_state_vector'] = '%15.6f   s' % time_of_first_state_vector
    time_of_second_state_vector = utc2sod(
        beijing2utc(GPSParam[1].find('TimeStamp').text))
    par['state_vector_interval'] = '%15.6f   s' % (
        time_of_second_state_vector - time_of_first_state_vector)
    state_vectors = [{i: float(node.find(j).text) for i, j in (('posX', 'xPosition'), ('posY', 'yPosition'), (
        'posZ', 'zPosition'), ('velX', 'xVelocity'), ('velY', 'yVelocity'), ('velZ', 'zVelocity'))} for node in
        GPSParam]

    keys = get_keys(True)
    res = ""
    res += ('Gamma Interferometric SAR Processor (ISP) - Image Parameter File\n\n')
    for k, v in zip(keys, par.values()):
        res += ('%s%s\n' % (k, v))
    for i, j in enumerate(state_vectors):
        res += ('state_vector_position_%d: %14.4f  %14.4f  %14.4f   m   m   m\n' % (
            i + 1, j['posX'], j['posY'], j['posZ']))
        res += ('state_vector_velocity_%d: %14.5f  %14.5f  %14.5f   m/s m/s m/s\n' % (
            i + 1, j['velX'], j['velY'], j['velZ']))
    return res


def th_cos2slc(cos_path):
    """
    将天绘二号卫星影像的数据文件.cos转换为gamma通用的数据文件slc格式

    :param cos_path: .cos文件的读取路径
    :return: slc
    """
    endian = 'little'  # TerraSAR-X的.cos是以大端序保存的，但天会二号是小端序。
    with open(cos_path, 'rb') as cos:
        s = 4
        Bytes_In_Burst = int.from_bytes(
            cos.read(s), byteorder=endian)  # BIB, 一个burst的总字节数
        Range_Sample_Relative_Index = int.from_bytes(
            cos.read(s), byteorder=endian)
        # RSRI, 距离向采样点索引值，给出了在一个虚拟栅格上每一行距离向上第一个采样点相对于图像中心参考点零多普勒时刻位置的索引值
        # RS, 距离向采样点数，该值在所有burst中保持一致
        Range_Samples = int.from_bytes(cos.read(s), byteorder=endian)
        Azimuth_Samples = int.from_bytes(
            cos.read(s), byteorder=endian)  # AS, 方位向采样点数
        Burst_Index = int.from_bytes(
            cos.read(s), byteorder=endian)  # BI, burst个数的索引值
        Rangeline_Total_Number_of_Bytes = int.from_bytes(
            cos.read(s), byteorder=endian)  # RTNB, 一行距离向的总字节数，包含注释信息
        Total_Number_of_Lines = int.from_bytes(cos.read(s), byteorder=endian)
        # TNL, 方位向上的长度，即行数，包含注释信息。该值只在文件第一行给出，方便文件读取和其他burst替换。因此整个文件大小可以通过RTNB×TNL得出
        # print('Bytes_In_Burst:', Bytes_In_Burst)
        # print('Range_Sample_Relative_Index:', Range_Sample_Relative_Index)
        # print('Range_Samples:', Range_Samples)
        # print('Azimuth_Samples:', Azimuth_Samples)
        # print('Burst_Index:', Burst_Index)
        # print('Rangeline_Total_Number_of_Bytes:',
        #       Rangeline_Total_Number_of_Bytes)
        # print('Total_Number_of_Lines:', Total_Number_of_Lines)
        # 很遗憾，上述数据中的BIB和RTNB是错误的，TNL也说不上对但暂且不管，并非读取错误而是数据本身错误。需要进行修正
        annotation_lines, annotation_columns, offset = 4 + 1, 2, -8
        # 天绘cos数据内的元数据是开头5行而不是开头4行，第5行还少了8字节，导致元数据的形状不是一个矩形
        Rangeline_Total_Number_of_Bytes = s * \
            (Range_Samples + annotation_columns)  # 真正的RTNB
        # 原始的RTNB = 2 * 真RTNB - 8，等于说除了Annotation（每行2字节）外长度增大至2倍，可能是设计者对short型理解有误
        # BIB = 原始的RTNB * TNL，因此BIB跟着RTNB出错
        # 真正的体积为 (TNL * 真RTNB) + (真RTNB - 8)，并不是一个矩形。
        cos.seek(annotation_lines * Rangeline_Total_Number_of_Bytes + offset)
        dtype = np.dtype(np.int16).newbyteorder(endian)
        res = np.float32(np.frombuffer(cos.read(Rangeline_Total_Number_of_Bytes * Azimuth_Samples), dtype=dtype).astype(
            np.dtype(np.int16).newbyteorder('big')).reshape(Azimuth_Samples, (Range_Samples + annotation_columns) * 2)[
            :, 4:])
        return res[:, ::2]+res[:, 1::2]*1j


def s1_bursts(xml_path):
    return int(etree.parse(xml_path).getroot().find("swathTiming").find("burstList").get("count"))


def s1_xml2par(xml_path, eof_path, burst):
    tree = etree.parse(xml_path)
    root = tree.getroot()
    imageAnnotation = root.find("imageAnnotation")
    imageInformation = imageAnnotation.find("imageInformation")
    slantRangeTime = imageInformation.find("slantRangeTime").text
    near_range_slc = float(slantRangeTime) * Const.c / 2
    azimuthTimeInterval = imageInformation.find("azimuthTimeInterval").text
    azimuth_line_time = float(azimuthTimeInterval)
    rangePixelSpacing = imageInformation.find("rangePixelSpacing")
    range_pixel_spacing = float(rangePixelSpacing.text)

    swathTiming = root.find("swathTiming")
    linesPerBurst = int(swathTiming.find("linesPerBurst").text)
    samplesPerBurst = int(swathTiming.find("samplesPerBurst").text)
    azimuthTime = swathTiming.find("burstList")[burst].find(
        "azimuthTime")  # azimuthTime
    par = {
        "line_offset": linesPerBurst * burst,
        "start_time": utc2sod(azimuthTime.text),
        "azimuth_line_time": azimuth_line_time,
        "near_range_slc": near_range_slc,
        "range_pixel_spacing": range_pixel_spacing,

        'range_samples': samplesPerBurst,
        'azimuth_lines': linesPerBurst,
        'image_format': "FCOMPLEX"
    }

    eof = etree.parse(eof_path).getroot()
    Validity_Start = eof.find("Earth_Explorer_Header").find(
        "Fixed_Header").find("Validity_Period").find("Validity_Start").text
    time_of_first_state_vector = utc2sod(Validity_Start.split("=")[1])
    time_of_first_state_vector -= 86400  # 开头第一个数据对应的应该是上一天的，未经过验证，可能出错
    state_vector_position, state_vector_velocity = [], []
    for osv in eof.find("Data_Block").find("List_of_OSVs"):
        state_vector_position.append(
            [float(osv.find("X").text), float(osv.find("Y").text), float(osv.find("Z").text)])
        state_vector_velocity.append(
            [float(osv.find("VX").text), float(osv.find("VY").text), float(osv.find("VZ").text)])
    par.update({
        "time_of_first_state_vector": time_of_first_state_vector,
        "state_vector_interval": 10.0,
        "state_vector_position": state_vector_position,
        "state_vector_velocity": state_vector_velocity
    })
    return par


def s1_tif2slc(tif_path):
    with TiffFile(tif_path) as img:
        return img.pages[0].asarray()
