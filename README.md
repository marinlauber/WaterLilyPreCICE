# WaterLily-PreCICE

WaterLily.jl bindings to connect to the PreCICE coupling library

## Usage

This package is a work in progress and is not yet available in the julia package registry. To use it, clone the repository and add the package to your julia environment:

```julia
] dev /path/to/WaterLilyPreCICE.jl
] instantiate
```

Then, you can use the package in your julia scripts:

```julia
using WaterLilyPreCICE
```

Running the different solvers

run `Julia`
```bash
julia --project=../../ Fluid.jl ../precice-config.xml
# preCICEJulia --project=/home/marin/Workspace/WaterLilyPreCICE/examples Fluid.jl precice-config.xml
```

run `G+Smo`
```bash
./Solid -c precice-config.xml
```

run `Calculix`
```bash
ccx_preCICE -i calculix -bg -precice-participant Calculix
```

## Examples

### Turek-Hron with G+Smo

...

![](assets/turek-hron.gif)


### Limped model interface