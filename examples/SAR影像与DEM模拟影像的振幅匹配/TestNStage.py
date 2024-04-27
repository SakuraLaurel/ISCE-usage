from isce.components.isceobj.StripmapProc.runRefineSecondaryTiming import fitOffsets
from isce.components.iscesys.Component.ProductManager import ProductManager
from isce.components.mroipac.ampcor.NStage import NStage
from isce.components import isceobj
import shelve

def estimateOffsetField(reference, secondary, azoffset=0, rgoffset=0):
    sim = isceobj.createSlcImage()
    sim.load(secondary+'.xml')
    sim.setAccessMode('READ')
    sim.createImage()

    sar = isceobj.createSlcImage()
    sar.load(reference + '.xml')
    sar.setAccessMode('READ')
    sar.createImage()

    objOffset = NStage(name='reference_offset1')
    objOffset.configure()
    objOffset.setAcrossGrossOffset(rgoffset)
    objOffset.setDownGrossOffset(azoffset)
    objOffset.setWindowSizeWidth(128)
    objOffset.setWindowSizeHeight(128)
    objOffset.setSearchWindowSizeWidth(40)
    objOffset.setSearchWindowSizeHeight(40)


    nAcross = 60
    nDown = 60

    if not objOffset.numberLocationAcross:
        objOffset.setNumberLocationAcross(nAcross)

    if not objOffset.numberLocationDown:
        objOffset.setNumberLocationDown(nDown)

    objOffset.setFirstPRF(1.0)
    objOffset.setSecondPRF(1.0)
    objOffset.setImageDataType1('complex')
    objOffset.setImageDataType2('complex')

    objOffset.nstage(sar, sim)

    sar.finalizeImage()
    sim.finalizeImage()

    result = objOffset.getOffsetField()
    return result

def TestNStage():
    referenceSlc = "reference_slc/reference.slc"
    secondarySlc = "coregisteredSlc/coarse_coreg.slc"
    # secondarySlc = "secondary_slc/secondary.slc"

    field = estimateOffsetField(referenceSlc, secondarySlc)

    outShelveFile = "test/misreg"
    odb = shelve.open(outShelveFile)
    odb['raw_field']  = field
    shifts, cull = fitOffsets(field,azazOrder=0,
            azrgOrder=0,
            rgazOrder=0,
            rgrgOrder=0,
            snr=5.0)
    odb['cull_field'] = cull
    odb['azpoly'] = shifts[0]
    odb['rgpoly'] = shifts[1]
    odb.close()     

    pm = ProductManager()
    pm.configure()
    pm.dumpProduct(shifts[0], outShelveFile + '_az.xml')
    pm.dumpProduct(shifts[1], outShelveFile + '_rg.xml')


if __name__ == "__main__":
    TestNStage()