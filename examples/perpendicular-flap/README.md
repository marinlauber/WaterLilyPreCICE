### Perpendicular Flap

To start the Fluid (`WaterLily`) participant

```bash
module load precice/3.1.1
preCICEJulia --project=../ Fluid.jl precice-config.xml
```

To start the Solid (`G+Smo`) participant

```bash
module load precice/3.1.1
make -C ~/Workspace/gismo/build/bin/test-2d-vertical-beam-implicit-solid
cp ~/Workspace/gismo/build/bin/test-2d-vertical-beam-implicit-solid Solid
./Solid -c precice-config.xml
```