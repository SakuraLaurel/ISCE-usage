from isce.components.iscesys.StdOEL.StdOELPy import create_writer
from isce.components.isceobj.Image.SlcImage import SlcImage
from isce.components.mroipac.ampcor.Ampcor import Ampcor
from isce.components.isceobj.Image.Image import Image
from isce.components.isceobj import createOffoutliers
from isce.applications.stripmapApp import Insar

def stripmapApp():
    insar = Insar(name = "stripmapApp",cmdline=["/home/sakura/Lys/test/src/stripmapApp.xml"])
    insar.configure()
    status = insar.run()
    raise SystemExit(status)

def main(test=False,windowSize=256,searchSize=16):
    sim_path = "test/sim.rdr"
    azoffset, rgoffset = 0, 0
    sim_type = "real"
    if test:
        sim_path = "secondary_slc/secondary.slc"
        azoffset, rgoffset = -27, 156
        sim_type = "mag"
    referenceSlc = "reference_slc/reference.slc"
    nAcross = 60
    nDown = 60
    azrgOrder=0
    azazOrder=0
    rgrgOrder=0
    rgazOrder=0

    sim = Image()
    sim.load(sim_path + '.xml')
    sim.setAccessMode('READ')
    sim.createImage()

    sar = SlcImage()
    sar.load(referenceSlc + '.xml')
    sar.setAccessMode('READ')
    sar.createImage()

    objOffset = Ampcor(name='')
    objOffset.configure()
    objOffset.setAcrossGrossOffset(rgoffset)
    objOffset.setDownGrossOffset(azoffset)
    objOffset.setWindowSizeWidth(windowSize)
    objOffset.setWindowSizeHeight(windowSize)
    objOffset.setSearchWindowSizeWidth(searchSize)
    objOffset.setSearchWindowSizeHeight(searchSize)
    margin = 2*objOffset.searchWindowSizeWidth + objOffset.windowSizeWidth

    nearestDistance = 21
    offAc = max(nearestDistance,-rgoffset)+margin
    offDn = max(nearestDistance,-azoffset)+margin

    lastAc = int( min(sar.getWidth(), sim.getWidth() - offAc) - margin)
    lastDn = int( min(sar.getLength(), sim.getLength() - offDn) - margin)

    objOffset.setFirstSampleAcross(offAc)
    objOffset.setLastSampleAcross(lastAc)
    objOffset.setFirstSampleDown(offDn)
    objOffset.setLastSampleDown(lastDn)
    objOffset.setNumberLocationAcross(nAcross)
    objOffset.setNumberLocationDown(nDown)
    objOffset.setFirstPRF(1.0)
    objOffset.setSecondPRF(1.0)
    objOffset.setImageDataType1('mag')
    objOffset.setImageDataType2(sim_type)

    objOffset.ampcor(sar, sim)

    sar.finalizeImage()
    sim.finalizeImage()

    field = objOffset.getOffsetField()

    aa, _ = field.getFitPolynomials(azimuthOrder=azazOrder, rangeOrder=azrgOrder, usenumpy=True)
    _, rr = field.getFitPolynomials(azimuthOrder=rgazOrder, rangeOrder=rgrgOrder, usenumpy=True)

    azshift = aa._coeffs[0][0]
    rgshift = rr._coeffs[0][0]
    print('Estimated az shift: ', azshift)
    print('Estimated rg shift: ', rgshift)

    import pickle
    with open("test/field60.pkl", "wb") as f:
        pickle.dump(field,f)


def filt(field, snr=5):
    for distance in [100,5,3,1]:
        objOff = createOffoutliers()
        objOff.wireInputPort(name='offsets', object=field)
        objOff.setSNRThreshold(snr)
        objOff.setDistance(distance)
        objOff.setStdWriter(create_writer("log","",True,filename='off.log'))
        objOff.offoutliers()
        field = objOff.getRefinedOffsetField()
        print('%d points left'%(len(field._offsets)))
        break

    aa, _ = field.getFitPolynomials(azimuthOrder=0, rangeOrder=0, usenumpy=True)
    _, rr = field.getFitPolynomials(azimuthOrder=0, rangeOrder=0, usenumpy=True)

    azshift = aa._coeffs[0][0]
    rgshift = rr._coeffs[0][0]
    print('Estimated az shift: ', azshift)
    print('Estimated rg shift: ', rgshift)

    # np.array(field.unpackOffsets())可以得到一个n x 5的数组，每一列依次是距离向采样位置、距离向偏移值、方位向采样位置、方位向偏移值、该次采样SNR。可以绘制出来，自己进行回归等操作。

# main(True)
# main(False)
if __name__ == "__main__":
    # main(False)
    import pickle
    with open("test/field.pkl","rb") as f:
        field = pickle.load(f)
        filt(field)