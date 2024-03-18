
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

def main(test=False,azrgOrder=0,azazOrder=0, rgrgOrder=0,rgazOrder=0,snr=5.0):
    sim_path = "test/sim.rdr"
    referenceSlc = "reference_slc/reference.slc"
    azoffset, rgoffset = 0, 0
    if test:
        sim_path = "test/sim.rdr"

    sim = Image()
    sim.load(sim_path + '.xml')
    sim.setAccessMode('READ')
    sim.createImage()

    sar = SlcImage()
    sar.load(referenceSlc + '.xml')
    sar.setAccessMode('READ')
    sar.createImage()

    width = sar.getWidth()
    length = sar.getLength()

    objOffset = Ampcor(name='reference_offset1')
    objOffset.configure()
    objOffset.setAcrossGrossOffset(rgoffset)
    objOffset.setDownGrossOffset(azoffset)
    objOffset.setWindowSizeWidth(128)
    objOffset.setWindowSizeHeight(128)
    objOffset.setSearchWindowSizeWidth(40)
    objOffset.setSearchWindowSizeHeight(40)
    margin = 2*objOffset.searchWindowSizeWidth + objOffset.windowSizeWidth

    nAcross = 60
    nDown = 60

    offAc = max(101,-rgoffset)+margin
    offDn = max(101,-azoffset)+margin

    lastAc = int( min(width, sim.getWidth() - offAc) - margin)
    lastDn = int( min(length, sim.getLength() - offDn) - margin)

    if not objOffset.firstSampleAcross:
        objOffset.setFirstSampleAcross(offAc)

    if not objOffset.lastSampleAcross:
        objOffset.setLastSampleAcross(lastAc)

    if not objOffset.firstSampleDown:
        objOffset.setFirstSampleDown(offDn)

    if not objOffset.lastSampleDown:
        objOffset.setLastSampleDown(lastDn)

    if not objOffset.numberLocationAcross:
        objOffset.setNumberLocationAcross(nAcross)

    if not objOffset.numberLocationDown:
        objOffset.setNumberLocationDown(nDown)

    objOffset.setFirstPRF(1.0)
    objOffset.setSecondPRF(1.0)
    objOffset.setImageDataType1('mag')
    objOffset.setImageDataType2('real')

    objOffset.ampcor(sar, sim)

    sar.finalizeImage()
    sim.finalizeImage()

    field = objOffset.getOffsetField()
    # for distance in [10,5,3,1]:
    #     inpts = len(field._offsets)
    #     print("DEBUG %%%%%%%%")
    #     print(inpts)
    #     print("DEBUG %%%%%%%%")
    #     objOff = createOffoutliers()
    #     objOff.wireInputPort(name='offsets', object=field)
    #     objOff.setSNRThreshold(snr)
    #     objOff.setDistance(distance)
    #     objOff.setStdWriter(create_writer("","",False))

    #     objOff.offoutliers()

    #     field = objOff.getRefinedOffsetField()
    #     outputs = len(field._offsets)

    #     print('%d points left'%(len(field._offsets)))

    aa, dummy = field.getFitPolynomials(azimuthOrder=azazOrder, rangeOrder=azrgOrder, usenumpy=True)
    dummy, rr = field.getFitPolynomials(azimuthOrder=rgazOrder, rangeOrder=rgrgOrder, usenumpy=True)

    azshift = aa._coeffs[0][0]
    rgshift = rr._coeffs[0][0]
    print('Estimated az shift: ', azshift)
    print('Estimated rg shift: ', rgshift)

main()