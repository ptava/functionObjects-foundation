## OpenFOAM function objects
---

This repository contains a collection of additional OpenFOAM foundation function objects.

---

### 1. refinementInfo
This function object writes any user-specified `uniformDimensionedScalarField` found in the registry at each `execute()` call.

Notes:
- created to write run-time into a file all info available associated with mesh refinement process;
- it needs modified `refiner` (see [modified dynamicFvMesh](https://github.com/ptava/fvMeshTopoChangers.git)) because is designed to access stored values in the registry such as `nTotToRefine`, `nTotRefined`, `lowerLimit`, `upperLimit` etc.

---

### 2. sasRefineIndicator
This function object generates a `volScalarField` that drives adaptive mesh refinement.

Notes:
- it is tailored for Scale-Adaptive Simulation (SAS) modelling, where mesh refinement should focus on regions of interest defined bu the turbulence length scales;

- requires length scale `Lvk` and its terms `C1` and `C2` (where $Lvk=max(C1,C2)$): available options to store these fields are in a modified `kOmegaSSTSAS` turbulence model (see [modified kOmegaSSTSAS](https://github.com/ptava/kOmegaSSTSAS.git)).

#### Available user options (refinement criteria):

Refinement indicator based on the **von Kármán length scale** ($L_{vk}$) and is
designed to serve as an adaptive refinement indicator for a Scale Adaptive
(SA) turbulente model.

Where:
- \f$ L_{vk} = \max(c_1, c_2) \f$
- \f$ c_1 \f$ is the physical turbulence length scale
- \f$ c_2 \f$ is the grid-related high wave number damper (linked to grid size)

The porpose is to adaptively refine mesh in regions where vortices are being
detected by a scale-resolving simulation, refining in vortices-associated
regions and coarsening outside them.

Options are provided to enphasise different regions of the flow/vortices:

1. **focusRegion == core**

Focuses refinement where the grid filter is most active, ideally in the
core of detected eddies.

Available transfer functions:

```
a) **markCoreConstant**:
```

Apply constant value to all cells where d = c2 - c1 > 0
\f[
    f(d) = \text{const}
\f]

```
b) **markCoreOddScaler**:
```

Apply an odd, monotonic, sign-preserving function, symmetric around
origin, normalised by its maximum value:

 <table align="center">
  <tr>
    <td align="center">
      <img src="imgs/tf_2.png" width="320" />
    </td>
    <td align="center">
      <img src="imgs/tf_3.png" width="320" />
    </td>
  </tr>
  <tr>
    <td align="center" colspan="2">
      <em>Transfer function representation with changing $w$ (left) and $\sigma$ (right)</em>
    </td>
  </tr>
</table>


\f[
  f(d) = \overline{d} \cdot \left( c - \exp\left( -\frac{\overline{d}^2}{2 \sigma^2} \right) \right), \quad c \geq 1.0
\f]

> [!NOTE]:
> - As \f$ |d| \to \infty \f$, the exponential vanishes, and \f$ f(d) \sim c \cdot d \f$
> - \f$ c = 1.0 \f$, function behaves cubically near 0 - useful for smoothly marking the interface region
> - \f$ \sigma \in [0.0, 1.0] \f$ controls the width of the transition region leading to steep decrease of \f$ f(d) \f$ near \f$ d^{max} \f$ (default: 0.0)
> - \f$ \overline{d} = d + w \cdot d^{max} \f$, with \f$ w \geq 0 \f$ - shifts the function to the left including in the refinement region also cells with small negative \f$ d \f$
> - f(d) is normalised by its maximum value to ensure consistent marking strength

<p align="center">
    <img src="imgs/os_1_0_18931.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 1$ and $w \eq 0</em>
</p>
<p align="center">
    <img src="imgs/os_1_0p3_22534.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 1$ and $w \eq 0.3</em>
</p>
<p align="center">
    <img src="imgs/os_1_0p5_24348.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 1$ and $w \eq 0.5</em>
</p>
<p align="center">
    <img src="imgs/os_1_1_28106.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 1$ and $w \eq 1</em>
</p>


```
c) **markCoreGaussSink**:
 ```

Gaussian-like function peaking where nLvk = Lvk / c2 = 1, and decays for
increasing values of nLvk:

 <table align="center">
  <tr>
    <td align="center">
      <img src="imgs/tf_0.png" width="320" />
    </td>
    <td align="center">
      <img src="imgs/tf_1.png" width="320" />
    </td>
  </tr>
  <tr>
    <td align="center" colspan="2">
      <em>Transfer function representation with changing $w$ (left) and $\sigma$ (right)</em>
    </td>
  </tr>
</table>

\f[
    f(nLvk) = \exp\left( -\frac{1}{2} \left( \frac{nLvk - 1}{\sigma} \right)^2 \right) + w \cdot (nLvk - 1)^2
\f]

<p align="center">
    <img src="imgs/gs_0p05_1_20037.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 0.05$ and $w \eq 1</em>
</p>
<p align="center">
    <img src="imgs/gs_0p5_1_24103.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 0.5$ and $w \eq 1</em>
</p>
<p align="center">
    <img src="imgs/gs_0p05_1_span.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 0.05$ and $w \eq 1</em>
</p>

<p align="center">
    <img src="imgs/gs_0p5_10_span.png" width="50%" height="50%">
    <em>Example with $\sigma \eq 0.05$ and $w \eq 10</em>
</p>


> [!NOTE]:
> - \f$ w \geq 0 \f$
> - \f$ w \f$ controls values outside the focus region
> - by construction, \f$ L_{vk} = max(c_1, c_2) \f$ therefore all cells where \f$ c_1 \leq c_2 \f$ (i.e., d >= 0) are going to have \f$ nLvk \leq 1 \f$

2. **focusRegion = "periphery"**

Focuses refinement in the periphery of vortex structures rather than
their core.

Available transfer functions:

```
d) **markPeripheryGaussSink**:
```

Gaussian-like function peaking at a user-defined reference value
of Lvk (e.g., a free-stream value). Same as 'c' but different nLvk
definition:

\f[
    f(nLvk) = w_1 \cdot \exp\left( -\frac{1}{2} \left( \frac{nLvk - 1}{\sigma} \right)^2 \right) - w_2 \cdot (nLvk - 1)^2
\f]


