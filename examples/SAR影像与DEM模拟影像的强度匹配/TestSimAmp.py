from isce.components.isceobj.Util.Simamplitude import Simamplitude
from isce.components.isceobj import createImage
from isce.components.iscesys.StdOEL.StdOELPy import create_writer

def simAmp():
    dem_path = "geometry/z.rdr.full.xml"
    sim_path = "test/sim.rdr"  # 注意不能取.amp后缀名，否则会被mdx当作是双通道数据
    
    hgtImage = createImage()
    hgtImage.load(dem_path)
    hgtImage.setAccessMode('read')
    hgtImage.createImage()

    simImage = createImage()
    simImage.setFilename(sim_path)
    simImage.dataType = 'FLOAT'
    simImage.setAccessMode('write')
    simImage.setWidth(hgtImage.getWidth())
    simImage.createImage()

    objShade = Simamplitude()
    stdWriter = create_writer("log","",False)
    objShade.setStdWriter(stdWriter)
    objShade.simamplitude(hgtImage, simImage, shade=3.0)  # shade=3.0是将模拟影像的强度增大到计算值的3倍

    simImage.renderHdr()
    hgtImage.finalizeImage()
    simImage.finalizeImage()

if __name__ == "__main__":
    simAmp()