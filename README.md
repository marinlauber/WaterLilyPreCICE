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

## Approach

The approach we have to use to couple a structural solver with `WaterLily.jl` using `preCICE` is different from what you would do with classical ALE methods. Since `WaterLily.jl` is an immersed-boundary solver, there is no interface mesh readily available in the fluid solver. We follow two approaches in this package; the first one is to read the mesh as triangulated surfaces and use the `preCICE` adapter to couple this mesh to whatever type of mesh is used in the structural solver. The second approach, which can only be used if we couple to `CalculiX` or `G+Smo` is to use `preCICE` to communicate that interface mesh directly from the structural solver to the fluid solver. The limitation to this is that only 

#### General coupling

For example, in the following figure, the left side represent the surface mesh required in `WaterLily.jl` and the right side the 3D mesh required in the structural solver. `WaterLily` uses this surface mesh to immerse the geometry and compute the forces on the faces of these elements. These forces are then passed to the structural solver (in this case [CalculiX](https://www.calculix.de/)). Here, we can use a nearest neighbor interpolation in the coupling since there is an exact correspondence between nodes/elements ion the two solvers.

 Surface mesh              |  3D mesh
:-------------------------:|:-------------------------:
![](assets/surface.png)    |  ![](assets/geometry.png)

This is the preferred approach if you structural problem is 3D (not a shell/membrane). 

> **_NOTE:_** 
Ths approach is also the one you want to be using to couple

#### Communicating mesh internally

If you structural problem consist of membrane/shells and you are using `CalculiC` or `G+Smo`m you can directly read the mesh information from the other participant.

<!-- > **_💡 Tip:_** -->
<!-- tiptip -->

## Examples

<!-- ### Rigid Sphere in a flow

...

![](assets/rigid-sphere.gif) -->


### Rubber Airfoil in a flow

This does not sound like a good idea, but it can actually help with gust mitigations!
![](assets/Airfoil_forces.png)

![](assets/rubber-airfoil.gif)


### Turek-Hron with G+Smo

...

![](assets/turek-hron.gif)


<!-- ### Lumped model interface -->

### Contributing

We always appreciate new contributions, so please submit a pull request with your changes and help us make this adapter better!

Of course, ideas, suggestions, and questions are welcome too! Please raise an issue to address any of these.

### Citing

