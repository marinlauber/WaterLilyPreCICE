import argparse
import ccx2paraview
import os

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
    os.system("rm -rf %s.12d %s.cvg %s.dat %s.sta %s.pvd datp" % (fn,fn,fn,fn,fn))
elif args.purge:
    print("Purging the files for %s" % filename)
    os.system("rm -rf %s.12d %s.cvg %s.dat %s.sta %s.pvd datp" % (fn,fn,fn,fn,fn))
    os.system("rm -rf %s.frd %s.inp %s.pvd %s.sta *.nam" % (fn,fn,fn,fn))
else:
    print("Post-processing file %s" % filename)
    c = ccx2paraview.Converter(filename, ["vtu"])
    c.run()
os.system("rm -rf WarnNodeMissMultiStage.nam spooles.out") 