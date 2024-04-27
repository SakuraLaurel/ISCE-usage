from os.path import getsize, basename, join
from tifffile import TiffFile  # type: ignore
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import re


class Slc(object):
    wavelength = 0.05546576
    dR = 2.329562114715323

    def __init__(self):
        self.referenceDir = "/home/sakura/Lys/data/S1A_IW_SLC__1SDV_20171104T101335_20171104T101402_019114_020569_D0D2.SAFE/measurement"
        self.secondaryDir = "/home/sakura/Lys/20181111"
        self.polarization = "vv"
        self.swath = 2
        self.burst = 1

    def display(self, attr="data", threshold=1):
        data: np.ndarray = getattr(self, attr)
        if data.dtype == np.complex64:
            data = np.abs(data)
        up = np.nanpercentile(data, 100 - threshold)
        data[data > up] = up
        down = np.nanpercentile(data, threshold)
        data[data < down] = down
        plt.imshow(data)
        plt.colorbar()
        plt.show()


class Filename(object):
    @classmethod
    def reference(cls, directory: str, swath: int, polarization: str) -> str:
        for i in glob(join(directory, "*")):
            if re.match(
                "[\s\S]*?iw%d[\s\S]*?%s[\s\S]*?" % (swath, polarization), basename(i)
            ):
                return i

    @classmethod
    def secondary(cls, directory: str, swath: int, burst: int) -> str:
        return join(directory, "fine_coreg/IW%d/burst_%02d.slc" % (swath, burst))

    @classmethod
    def rgOffset(cls, directory: str, swath: int, burst: int) -> str:
        return join(directory, "fine_offsets/IW%d/range_%02d.off" % (swath, burst))

    @classmethod
    def interferogram(cls, directory: str, swath: int, burst: int) -> str:
        return join(
            directory, "fine_interferogram/IW%d/burst_%02d.int" % (swath, burst)
        )

    @classmethod
    def multilook(cls, directory: str, swath: int, burst: int) -> str:
        return join(
            directory,
            "fine_interferogram/IW%d/burst_%02d.7alks_19rlks.int" % (swath, burst),
        )

    @classmethod
    def correlation(cls, directory: str, swath: int, burst: int) -> str:
        return join(
            directory, "fine_interferogram/IW%d/burst_%02d.cor" % (swath, burst)
        )


class Reference(Slc):
    @property
    def swath(self) -> int:
        return self._swath

    @swath.setter
    def swath(self, value: int):
        self._swath = value
        path = Filename.reference(self.referenceDir, value, self.polarization)
        with TiffFile(path) as img:
            self._data = img.pages[0].asarray()
            self._shape = self._data.shape
            self.length = int(self._shape[0] / 9)
            self.width = self._shape[1]
            self.shape = (self.length, self.width)  # (1506, 25441)

    @property
    def data(self):
        return self._data[(self.burst - 1) * self.length : self.burst * self.length]


class Secondary(Slc):
    def __init__(self, ref: Reference):
        super().__init__()
        self.ref = ref.data
        for i in ("length", "width", "shape"):
            setattr(self, i, getattr(ref, i))
        path = Filename.secondary(self.secondaryDir, self.swath, self.burst)
        data = self._fromfile(path, np.complex64, self.shape)
        path = Filename.rgOffset(self.secondaryDir, self.swath, self.burst)
        rgOffset = self._fromfile(path, np.float32, self.shape)
        phase = np.exp(1j * 4 * np.pi * rgOffset * self.dR / self.wavelength)
        self.data = data * phase

    @property
    def interferogram(self):
        path = Filename.interferogram(self.secondaryDir, self.swath, self.burst)
        return self._fromfile(path, np.complex64, self.shape)

    @property
    def correlation(self):
        path = Filename.correlation(self.secondaryDir, self.swath, self.burst)
        data = np.fromfile(path, np.float32).reshape((self.length * 2, self.width))
        mag, cor = data[::2], data[1::2]
        return cor

    def checkInterferogram(self):
        interferogram = self.ref * np.conj(self.data)
        return self._check(self.interferogram, interferogram)

    def checkMultiLook(self) -> np.ndarray:
        path = Filename.multilook(self.secondaryDir, self.swath, self.burst)
        match = re.search(r"(\d+?)alks_(\d+?)rlks", path)
        az = int(match.group(1))
        rg = int(match.group(2))
        length = int(self.length / az)
        width = int(self.width / rg)
        shape = (length, width)
        multilook = np.zeros(shape, dtype=np.complex64)
        interferogram = self.interferogram
        for i in range(az):
            for j in range(rg):
                multilook += interferogram[i : az * length : az, j : rg * width : rg]
        multilook /= az * rg
        return self._check(self._fromfile(path, np.complex64, shape), multilook)

    def checkCorrelation(self) -> np.ndarray:
        return self.correlation - self._calCorrelation()

    def checkMagnitude(self) -> np.ndarray:
        path = Filename.correlation(self.secondaryDir, self.swath, self.burst)
        data = np.fromfile(path, np.float32).reshape((self.length * 2, self.width))
        mag = data[::2]
        return np.sqrt(np.abs(self.data) * np.abs(self.ref)) - mag

    def _fromfile(
        self, path: str, dtype: np.dtype, shape: tuple[int, int]
    ) -> np.ndarray:
        return np.fromfile(path, dtype).reshape(shape)

    def _check(self, i, j):
        magnitude = np.abs(j - i) / np.abs(i)
        angle = np.angle(i * np.conj(j))
        return magnitude - angle
        # return magnitude

    def _calCorrelation(self, windowSize=5) -> np.ndarray:
        half = int(windowSize // 2)
        w = np.array(
            [1 - np.abs(2 * (i - half) / (windowSize + 1)) for i in range(windowSize)]
        )
        multiply = self.ref * self.data.conj()
        ref2 = np.abs(self.ref) ** 2
        sec2 = np.abs(self.data) ** 2
        p1 = np.zeros(self.shape, np.complex64)
        p2 = np.zeros(self.shape, np.float32)
        p3 = np.zeros(self.shape, np.float32)
        res = np.zeros(self.shape, np.float32)
        resSlice = slice(half, 1 + half - windowSize)
        for i in range(windowSize):
            iSlice = slice(i, i + self.length + 1 - windowSize)
            for j in range(windowSize):
                jSlice = slice(j, j + self.width + 1 - windowSize)
                wt = w[i] * w[j]
                p1[resSlice, resSlice] += multiply[iSlice, jSlice] * wt
                p2[resSlice, resSlice] += ref2[iSlice, jSlice] * wt
                p3[resSlice, resSlice] += sec2[iSlice, jSlice] * wt
        res[resSlice, resSlice] = (
            np.abs(p1[resSlice, resSlice])
            / (p2[resSlice, resSlice] * p3[resSlice, resSlice]) ** 0.5
        )
        return res


if __name__ == "__main__":
    ref = Reference()
    sec = Secondary(ref)
    sec.temp = sec.checkMagnitude()[200:300, 3300:3400]
    sec.display("temp")
