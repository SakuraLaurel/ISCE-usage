from isce.components.isceobj.StripmapProc.runRefineSecondaryTiming import estimateOffsetField, fitOffsets
from isce.components.iscesys.Component.ProductManager import ProductManager
import shelve

def runRefineSecondaryTiming(beforeResample=True):
    secondarySlc = "coregisteredSlc/coarse_coreg.slc"
    rgoffset = 0
    if beforeResample:
        secondarySlc = "secondary_slc/secondary.slc"
        rgoffset = 150
    referenceSlc = "reference_slc/reference.slc"
    
    field = estimateOffsetField(referenceSlc, secondarySlc,0,rgoffset)

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
    runRefineSecondaryTiming(beforeResample=True)