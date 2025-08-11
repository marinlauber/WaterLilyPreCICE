# get G+Smo
git clone git@github.com:gismo/gismo
cd gismo
git checkout stable

# add optional modules
cd optional

# clone gsStructuralAnalysis
git clone git@github.com:gismo/gsStructuralAnalysis

# clone gsKLShell
git clone git@github.com:gismo/gsKLShell

# clone gsPreCICE and gsElasicity
git clone git@github.com:gismo/gsPreCICE
git clone git@github.com:gismo/gsElasticity

# configure the make files
cd ../
mkdir build && cd build
cmake -DGISMO_WITH_OPENMP=ON -DGISMO_WITH_MPI=ON -DGISMO_BUILD_EXAMPLES=OFF -DGISMO_OPTIONAL="gsSpectra;gsElasticity;gsKLShell;gsStructuralAnalysis;gsPreCICE" ../
make -j 8