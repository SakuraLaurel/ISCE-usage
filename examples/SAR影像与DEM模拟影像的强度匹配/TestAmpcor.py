from isce.components.isceobj.StripmapProc.runRefineSecondaryTiming import estimateOffsetField, fitOffsets
from isce.components.iscesys.Component.ProductManager import ProductManager
import shelve

def runRefineSecondaryTiming():
    pm = ProductManager()
    pm.configure()
    # secondarySlc = "coregisteredSlc/coarse_coreg.slc"
    secondarySlc = "secondary_slc/secondary.slc"
    field = estimateOffsetField("reference_slc/reference.slc", secondarySlc)

    rgratio = 1
    azratio = 1

    print ('*************************************')
    print ('rgratio, azratio: ', rgratio, azratio)
    print ('*************************************')

    outShelveFile = "test/misreg"
    odb = shelve.open(outShelveFile)
    odb['raw_field']  = field
    shifts, cull = fitOffsets(field,azazOrder=0,
            azrgOrder=0,
            rgazOrder=0,
            rgrgOrder=0,
            snr=5.0)
    odb['cull_field'] = cull

    ####Scale by ratio
    for row in shifts[0]._coeffs:
        for ind, val in  enumerate(row):
            row[ind] = val * azratio

    for row in shifts[1]._coeffs:
        for ind, val in enumerate(row):
            row[ind] = val * rgratio


    odb['azpoly'] = shifts[0]
    odb['rgpoly'] = shifts[1]
    odb.close()     

    pm.dumpProduct(shifts[0], outShelveFile + '_az.xml')
    pm.dumpProduct(shifts[1], outShelveFile + '_rg.xml')


if __name__ == "__main__":
    runRefineSecondaryTiming()