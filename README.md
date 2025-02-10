# WaterLily-PreCICE

[WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl) bindings to connect to the [PreCICE](https://precice.org) coupling library

## Prerequisites

Before running any coupled simulation using `WaterLily.jl` and any other package with an excisting preCICE adapter, you will need a working version of preCICE, see [here](https://precice.org/installation-overview.html) for a detailed installation guide.

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

### Rigid Sphere in a flow

...

![](assets/rigid-sphere.gif)

### Turek-Hron with G+Smo

...

![](assets/turek-hron.gif)


<!-- ### Lumped model interface -->

### Contributing

We always appreciate new contributions, so please submit a pull request with your changes and help us make this adapter better!

Of course, ideas, suggestions, and questions are welcome too! Please raise an issue to address any of these.

### Citing

