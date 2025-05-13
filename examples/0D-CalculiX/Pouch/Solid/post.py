import argparse
import ccx2paraview
import os

def _clean(fname="calculix", exts=["frd","12d","cvg","dat","sta","pvd"]):
    for ex in exts:
        os.system("rm -rf %s.%s" % (fn,ex))
    return None

# get the legnth scale if provided
parser = argparse.ArgumentParser(description="Post-process the Calculix.frd file")
parser.add_argument("-F","--file", help="The file name to pst-process")
parser.add_argument("-C","--clean", help="Clean the useless files",action="store_true")
parser.add_argument("-P","--purge", help="Purge all the results and pvd files",action="store_true")
args = parser.parse_args()

# set default variables
filename = "calculix.frd" if(not args.file) else args.file

# run the postprocessing
fn = filename.split(".")[0]
if args.clean:
    print("Cleaning the files for %s" % filename)
    _clean(filename)
elif args.purge:
    print("Purging the files for %s" % filename)
    _clean(fn); _clean(filename,["inp","nam"])
else:
    print("Post-processing file %s" % filename)
    c = ccx2paraview.Converter(filename, ["vtu"])
    c.run()
os.system("rm -rf WarnNodeMissMultiStage.nam spooles.out") 