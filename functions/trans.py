from datetime import datetime, timedelta
import numpy as np
import re


class Const(object):
    """
    常数
    """
    a = 6378137.0000  # 地球半长轴
    b = 6356752.3141  # 地球半短轴
    c = 299792458  # 光速，不是地球焦半径
    e2 = 1 - (b / a) ** 2  # 偏心率的平方
    ep2 = (a / b) ** 2 - 1  # 第二偏心率的平方
    f = 1 - b / a  # 扁率
    UTMScaleFactor = 0.9996  # UTM坐标系中央经线的比例因子
    n = (a - b) / (a + b)  # 我不知道这叫啥，但在坐标转换的时候要用


def degree2rad(degree):
    """
    角度制转弧度制

    :param degree: 角度值
    :return: 弧度值
    """
    return degree * np.pi / 180


def rad2degree(rad):
    """
    弧度制转角度制

    :param rad: 弧度值
    :return: 角度值
    """
    return rad * 180 / np.pi


def utc2sod(UTCtime):
    """
    UTC时间转换为日积秒
    """
    date = re.split('[-T:Z]', UTCtime)[:6]
    return int(date[3]) * 3600 + int(date[4]) * 60 + float(date[5])


def beijing2utc(time, format='%Y-%m-%dT%H', mark=':'):
    """
    天绘影像为北京时间，该函数将北京时间转换为UTC时间

    :param time: 北京时间，YYYY-MM-DDTHH:分钟:秒，冒号之后的无需修改
    :param format: 照着天绘数据格式给出的解析式
    :param mark: 冒号之前的需要修改，之后的照抄就行
    :return: UTC时间
    """
    t = time.split(mark)  # 分割，只需修改t[0]，即年月日时
    t[0] = (datetime.strptime(t[0], format) -
            timedelta(hours=8)).strftime(format)
    return mark.join(t)  # 补上分钟和秒


def cos_vectors(v1, v2, epsilon=1e-8):
    """
    计算两个向量夹角的余弦

    :param v1: 向量1
    :param v2: 向量2
    :return: cos<v1, v2>
    """
    if v1.ndim == 1:
        shape = (1,)
    else:
        shape = v1.shape[1:]
    v3, v4 = np.reshape(v1, (3, -1)), np.reshape(v2, (3, -1))
    n1, n2 = np.linalg.norm(v3, axis=0), np.linalg.norm(v4, axis=0)
    if np.any(n1 < epsilon) or np.any(n2 < epsilon):
        raise ValueError("零向量的夹角无法计算")
    res = np.sum(np.multiply(v3, v4), axis=0)/(n1 * n2)
    return res.reshape(shape)


def doy(time, format='%Y%m%d'):
    return (datetime.strptime(time, format) - datetime.strptime(time[:4] + "0101",format)).days + 1