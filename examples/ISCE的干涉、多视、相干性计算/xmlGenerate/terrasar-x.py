from os.path import basename
from lxml import etree
from glob import glob
from os import system
import re


def getDate(filename):
    return re.findall(r"_\d{8}T", basename(filename))[0][1:-1]


def filePath(filename: str):
    if filename == basename(filename):
        f = filename.split(".")[0]
        return "/home/sakura/laurel/DATA/TerraSAR-X/%s/TSX-1.SAR.L1B/%s/%s" % (
            f,
            f,
            filename,
        )
    else:
        return filename


def generate(reference, secondary):
    topsApp = etree.Element("insarApp")
    component = etree.SubElement(topsApp, "component", {"name": "insar"})
    property = etree.SubElement(component, "property", {"name": "sensor name"})
    property.text = "TERRASARX"
    component2 = etree.SubElement(component, "component", {"name": "reference"})
    property2 = etree.SubElement(component2, "property", {"name": "XML"})
    property2.text = filePath(reference)
    property3 = etree.SubElement(component2, "property", {"name": "OUTPUT"})
    property3.text = getDate(reference)

    component2 = etree.SubElement(component, "component", {"name": "secondary"})
    property2 = etree.SubElement(component2, "property", {"name": "XML"})
    property2.text = filePath(secondary)
    property3 = etree.SubElement(component2, "property", {"name": "OUTPUT"})
    property3.text = getDate(secondary)

    property = etree.SubElement(component, "property", {"name": "demfilename"})
    property.text = "/home/sakura/lys/tx/demLat_N39_N41_Lon_E116_E117.dem.wgs84"
    return etree.tostring(topsApp, encoding="utf-8", pretty_print=True).decode("utf-8")


# for i in glob("/home/sakura/laurel/DATA/TerraSAR-X/*/TSX-1.SAR.L1B/T*"):
#     source = re.match("(/home[\s\S]*?)/TSX-1.SAR.L1B", i).group(1)
#     target = source.replace(basename(source), basename(i))
#     system('mv "%s" "%s"' % (source, target))

reference = "TDX1_SAR__SSC______SM_S_SRA_20171114T100443_20171114T100453.xml"
for secondary in glob("/home/sakura/laurel/DATA/TerraSAR-X/*/TSX-1.SAR.L1B/T*/T*"):
    if not re.findall(r"20171114", secondary):
        res = generate(reference, secondary)
        with open("./sh/" + getDate(secondary) + ".xml", "w") as f:
            f.write(res)
