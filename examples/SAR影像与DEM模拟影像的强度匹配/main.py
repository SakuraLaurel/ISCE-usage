from isce.applications.stripmapApp import Insar

def main():
    insar = Insar(name = "stripmapApp",cmdline=["/home/sakura/Lys/test/src/stripmapApp.xml"])
    insar.configure()
    status = insar.run()
    raise SystemExit(status)

main()