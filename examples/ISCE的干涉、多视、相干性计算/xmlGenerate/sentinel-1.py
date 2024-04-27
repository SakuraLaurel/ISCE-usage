from os.path import basename
from glob import glob
import re
from lxml import etree


def getDate(filename):
    return re.findall(r"_\d{8}T", basename(filename))[0][1:-1]


def filePath(filename):
    if filename == basename(filename):
        return "/home/sakura/Lys/data/%s" % filename
    else:
        return filename


def auxDir():
    return "/home/sakura/Lys/S1A_AUX_CAL_V20171017T080000_G20201215T124601.SAFE"


def generate(reference, secondary):
    topsApp = etree.Element("topsApp")
    component = etree.SubElement(topsApp, "component", {"name": "topsApp"})
    property = etree.SubElement(component, "property", {"name": "sensor name"})
    property.text = "SENTINEL1"
    component2 = etree.SubElement(component, "component", {"name": "reference"})
    property2 = etree.SubElement(component2, "property", {"name": "safe"})
    value = etree.SubElement(property2, "value")
    value.text = filePath(reference)
    property3 = etree.SubElement(component2, "property", {"name": "orbit directory"})
    value2 = etree.SubElement(property3, "value")
    # value2.text = matchOrbit(reference)
    value2.text = "/home/sakura/Lys/orbits/"
    property4 = etree.SubElement(component2, "property", {"name": "output directory"})
    value3 = etree.SubElement(property4, "value")
    value3.text = "reference"
    property5 = etree.SubElement(
        component2, "property", {"name": "auxiliary data directory"}
    )
    value4 = etree.SubElement(property5, "value")
    value4.text = auxDir()
    # property6 = etree.SubElement(component2, "property", {"name": "swath number"})
    # value5 = etree.SubElement(property6, "value")
    # value5.text = "2"

    component2 = etree.SubElement(component, "component", {"name": "secondary"})
    property2 = etree.SubElement(component2, "property", {"name": "safe"})
    value = etree.SubElement(property2, "value")
    value.text = filePath(secondary)
    property3 = etree.SubElement(component2, "property", {"name": "orbit directory"})
    value2 = etree.SubElement(property3, "value")
    # value2.text = matchOrbit(reference)
    value2.text = "/home/sakura/Lys/orbits/"
    property4 = etree.SubElement(component2, "property", {"name": "output directory"})
    value3 = etree.SubElement(property4, "value")
    value3.text = "secondary"
    property5 = etree.SubElement(
        component2, "property", {"name": "auxiliary data directory"}
    )
    value4 = etree.SubElement(property5, "value")
    value4.text = auxDir()
    # property6 = etree.SubElement(component2, "property", {"name": "swath number"})
    # value5 = etree.SubElement(property6, "value")
    # value5.text = "2"

    property = etree.SubElement(component, "property", {"name": "demfilename"})
    property.text = "/home/sakura/Lys/dem/demLat_N38_N41_Lon_E114_E118.dem.wgs84"
    return etree.tostring(topsApp, encoding="utf-8", pretty_print=True).decode("utf-8")


if __name__ == "__main__":
    reference = (
        "S1A_IW_SLC__1SDV_20171104T101335_20171104T101402_019114_020569_D0D2.SAFE"
    )
    for secondary in glob("/home/sakura/Lys/data/*"):
        if not re.findall(r"20171104", secondary):
            res = generate(reference, secondary)
            with open("./sh/" + getDate(secondary) + ".xml", "w") as f:
                f.write(res)

# conda activate isce && cd /home/sakura/gran && topsApp.py /home/sakura/laurel/closure/sh/20171116.xml
