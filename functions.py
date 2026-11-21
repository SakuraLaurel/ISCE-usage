from glob import glob
from functools import reduce
from os.path import join, basename
from pickle import dump, load
from multiprocessing import Pool

def lmap(function, iterable):
    return list(map(function, iterable))

def lreduce(function, iterable):
    if len(iterable) <= 1:
        return list(iterable)
    else:
        return list(reduce(function, iterable))

def flip(func, *args):
    return func(*args[::-1])

def pmap(func, iter):
    p = Pool()
    return p.map(func, iter)

def find_files(filename):
    path = "/home/sakura/p6"
    dirs = ("/home/sakura/p1",
            "/home/sakura/p5",
            path,
            join(path, "aux/dem"),
            join(path, "aux/landcover"))
    res = lmap(lambda dir: glob(join(dir, filename)), dirs)
    return sorted(lreduce(lambda x, y: x + y, res), key=lambda x: basename(x))

lcs = {
    "forest": "森林",
    "grassland": "草地",
    "impervious": "城市",
    "cropland": "耕地"
}
region = {
    "forest": (40.58, 116.43, 40.51, 116.5),
    "grassland": (40.07, 115.7, 40, 115.77),
    "impervious": (40, 116.3,39.93,116.37),
    "cropland": (39.3, 116.63,39.23,116.7)
}
swath = {
    "forest": 3,
    "grassland": 2,
    "impervious": 3,
    "cropland": 3
}
pkl_dir = "data"


def dump_pkls(lc, date, data):
    with open(join(pkl_dir, f"{lc}/{date}.pkl"), "wb") as f:
        dump(data, f)

def load_pkls(lc, date):
    with open(join(pkl_dir, f"{lc}/{date}.pkl"), "rb") as f:
        return load(f)


# 这个东西是循环色带，很有趣，但是在本项目中不会用到
# from matplotlib.colors import LinearSegmentedColormap
# N_total = 256
# colors = [(1.0, 0.0, 0.0),  # 白色
#           (0.0, 0.4, 1.0),  # 蓝色
#           (1.0, 1.0, 1.0),  # 深红色
#           (0.0, 0.4, 1.0),  # 再次蓝色，为了平滑过渡
#           (1.0, 0.0, 0.0)]  # 白色
# cmap = LinearSegmentedColormap.from_list("housyou's cmap", colors, N=N_total )