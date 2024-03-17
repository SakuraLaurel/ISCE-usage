from isce.components.isceobj.StripmapProc.runDenseOffsets import estimateOffsetField

def runDenseOffsets(denseWindowWidth=64, denseWindowHeight=64,denseSearchWidth=20,denseSearchHeight=20,denseSkipWidth=32,denseSkipHeight=32):
    referenceSlc = "reference_slc/reference.slc"
    # secondarySlc = "secondary_slc/secondary.slc"
    secondarySlc = "coregisteredSlc/coarse_coreg.slc"
    denseOffsetFilename = "test/res"

    field = estimateOffsetField(referenceSlc, secondarySlc, denseOffsetFilename,
                                ww = denseWindowWidth,
                                wh = denseWindowHeight,
                                sw = denseSearchWidth,
                                shh = denseSearchHeight,
                                kw = denseSkipWidth,
                                kh = denseSkipHeight)

    print(field[0], field[1])

runDenseOffsets()