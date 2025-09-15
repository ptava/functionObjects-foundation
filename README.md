# OpenFOAM function objects
---

This repository contains a collection of additional OpenFOAM foundation function objects.

### 1. refinementInfo
This function object writes any user-specified `uniformDimensionedScalarField` found in the registry at each `execute()` call.

Notes:
- created to write run-time into a file all info available associated with mesh refinement process;
- it needs modified `refiner` (see [modified dynamicFvMesh](https://github.com/ptava/fvMeshTopoChangers.git)) because is designed to access stored values in the registry such as `nTotToRefine`, `nTotRefined`, `lowerLimit`, and `upperLimit`.

### 2) sasRefineIndicator
This function object generates a `volScalarField` that drives adaptive mesh refinement.

Notes:
- it is tailored for Scale-Adaptive Simulation (SAS) modelling, where mesh refinement should focus on regions of interest defined bu the turbulence length scales;

- this implementation requires a modified `kOmegaSSTSAS` turbulence model (see [modified kOmegaSSTSAS](https://github.com/ptava/kOmegaSSTSAS.git)) to make available the length scale `Lvk` its terms `C1` and `C2`, where $Lvk=max(C1,C2)$;

- Available user options/criteria:
    * `core`: it marks for refinement those cells where the high wave number damper term is active (i.e., for kOmegaSSTSAS model, where $C2 >= C1$). Under testing the way we can mark those cells, with constant or Gaussian function;
    * `periphery`: it marks for refinement those cells at the boundaries of the region where unsteadiness is detected. Under development;
    * `combined`: use both above-mentioned criteria to mark cells for refinement (not yet implemented).
