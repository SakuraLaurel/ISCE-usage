from os.path import basename
from lxml import etree
from glob import glob
from os import system
import re


def getDate(filename):
    return "20" + re.findall(r"-\d{6}", basename(filename))[0][1:]


def filePath(filename: str):
    if filename != basename(filename):
        return filename
    else:
        return glob("/home/sakura/laurel/DATA/ALOS-2/*/%s" % filename)[0]


def leaderfilePath(imagefile: str):
    imagefile = filePath(imagefile)
    return imagefile.replace("IMG-HH-", "LED-")


def generate(reference, secondary):
    topsApp = etree.Element("insarApp")
    component = etree.SubElement(topsApp, "component", {"name": "insar"})
    property = etree.SubElement(component, "property", {"name": "sensor name"})
    property.text = "ALOS2"
    component2 = etree.SubElement(component, "component", {"name": "reference"})
    property2 = etree.SubElement(component2, "property", {"name": "IMAGEFILE"})
    property2.text = filePath(reference)
    property3 = etree.SubElement(component2, "property", {"name": "LEADERFILE"})
    property3.text = leaderfilePath(reference)
    property4 = etree.SubElement(component2, "property", {"name": "OUTPUT"})
    property4.text = getDate(reference)

    component2 = etree.SubElement(component, "component", {"name": "secondary"})
    property2 = etree.SubElement(component2, "property", {"name": "IMAGEFILE"})
    property2.text = filePath(secondary)
    property3 = etree.SubElement(component2, "property", {"name": "LEADERFILE"})
    property3.text = leaderfilePath(secondary)
    property4 = etree.SubElement(component2, "property", {"name": "OUTPUT"})
    property4.text = getDate(secondary)

    property = etree.SubElement(component, "property", {"name": "demfilename"})
    property.text = "/home/sakura/lys/alos2/demLat_N39_N41_Lon_E116_E118.dem.wgs84"
    return etree.tostring(topsApp, encoding="utf-8", pretty_print=True).decode("utf-8")


reference = "IMG-HH-ALOS2293430790-191029-FBDR1.1__A"
for secondary in glob("/home/sakura/laurel/DATA/ALOS-2/*/IMG-HH*"):
    if not re.findall(r"191029", secondary):
        res = generate(reference, secondary)
        with open("./sh/" + getDate(secondary) + ".xml", "w") as f:
            f.write(res)
