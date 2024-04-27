from os import mkdir, system, chdir, cpu_count
from os.path import basename, exists, join
from multiprocessing import Pool
from datetime import datetime, timedelta
from glob import glob
import numpy as np
import re


def downEOF():
    ans = []
    with open("temp.txt", "r") as f:
        source = f.read()
        for i in glob("/home/sakura/Lys/data/*"):
            d = (
                r"S1A_OPER_AUX_POEORB_OPOD_[0-9_T]*?V"
                + theDayBefore(getDate(i))
                + r"T[0-9_T]*?.EOF"
            )
            ans.append(re.findall(d, source)[-1])
    with open("res.txt", "w") as f:
        f.write(
            "\n".join(["https://s1qc.asf.alaska.edu/aux_poeorb/%s" % j for j in ans])
        )


def theDayBefore(sDate):
    return (datetime.fromisoformat(sDate) - timedelta(1)).strftime("%Y%m%d")


def getDate(filename):
    return re.findall(r"_\d{8}T", basename(filename))[0][1:-1]


def matchOrbit(filename):
    for i in glob("/home/sakura/Lys/orbits/*.EOF"):
        if re.match(
            r"S1A_OPER_AUX_POEORB_OPOD_[0-9_T]*?V%sT[0-9_T]*?.EOF"
            % theDayBefore(getDate(filename)),
            basename(i),
        ):
            return i


def unzip():
    """
    多线程命令可能要整个复制到if __name__ == "__main__"下执行
    """
    chdir("/home/sakura/Lys/")
    n = cpu_count()
    files = sorted(glob("/home/sakura/Lys/Downloads/*.zip"))
    each = int(np.ceil(len(files) / n))
    p = Pool(n)
    for i in range(n):
        p.apply_async(
            lambda j: [system("unzip %s" % k) for k in j],
            args=(files[each * i : each * (i + 1)],),
        )
    p.close()
    p.join()


def sentinel():
    for i in sorted(glob("/home/sakura/laurel/closure/sh/*.xml"))[54:]:
        filename = basename(i).split(".")[0]
        path = join("/home/sakura/Lys/result", filename)
        if not exists(path):
            mkdir(path)
        system("cd %s && topsApp.py %s" % (path, i))


def terrasarx():
    for i in [1, 2, 16, 17]:
        j = sorted(glob("/home/sakura/laurel/closure/sh/*.xml"))[i]
        filename = basename(j).split(".")[0]
        path = join("/home/sakura/lys/tx", filename)
        if not exists(path):
            mkdir(path)
        system("cd %s && stripmapApp.py %s" % (path, j))


def alos2():
    for i in [0, 2, 3]:
        j = sorted(glob("/home/sakura/laurel/closure/alos2/*.xml"))[i]
        filename = basename(j).split(".")[0]
        path = join("/home/sakura/lys/alos2", filename)
        if not exists(path):
            mkdir(path)
        system("cd %s && stripmapApp.py %s" % (path, j))


if __name__ == "__main__":
    # sentinel()
    # terrasarx()
    # alos2()
    pass
