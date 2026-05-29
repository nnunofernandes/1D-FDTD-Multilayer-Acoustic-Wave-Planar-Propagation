# 1D FDTD Planar Acoustic Wave Propagation
[![DOI](https://zenodo.org/badge/1248269950.svg)](https://doi.org/10.5281/zenodo.20445079)

## Overview
This repository contains a Python implementation of a one-dimensional finite-difference time-domain (FDTD) model for simulating planar acoustic wave propagation through layered media. The code was developed to support the numerical study of high-frequency ultrasound propagation through biological tissues and transducer-coupled material systems.

The model is intended for exploratory simulation, parameter screening, and preliminary estimation of acoustic transmission through tissue-like layers. It can be adapted to different propagation paths by modifying the material stack, spatial discretization, excitation signal, and boundary conditions.

## Repository purpose

The script simulates the propagation of an acoustic pressure wave along a one-dimensional spatial domain. Each region of the domain can be assigned different acoustic properties, allowing the user to represent layered systems such as:

- piezoelectric transducer interface layers;
- coupling media;
- skin;
- adipose tissue;
- muscle;
- tendon or ligament;
- cartilage;
- cortical or trabecular bone;
- simplified tissue-mimicking phantoms.

<img width="3879" height="1661" alt="KneeFDTD" src="https://github.com/user-attachments/assets/635d8366-2d7a-4f44-97db-8946e95aff97" />

The model is particularly useful when the main objective is to estimate how much acoustic energy may be transmitted, reflected, delayed, or attenuated as a function of tissue composition and propagation distance.

## Scientific background

The implemented method follows the finite-difference time-domain approach for acoustic wave propagation. In this formulation, the continuous acoustic wave equation is discretized in space and time, allowing the pressure field to be updated iteratively.

<img width="4082" height="1744" alt="FDTD" src="https://github.com/user-attachments/assets/8f20a8b2-ed87-4cbc-8967-07666c6c5fed" />


In a simplified linear acoustic formulation, the pressure and velocity field may be represented by the system of couple equations:

```math
    c^{n+\frac{1}{2}}_i=
    c_i^{n-\frac{1}{2}}-
    (p_{i+1}^n -p_i^n)
    \cdot \frac{\Delta t}{\rho _i \cdot \Delta x} ,
```

```math
    p^{n+\frac{1}{2}}_i=
    c_i^n-
    (c_{i+1}^{n+\frac{1}{2}} -c_i^{n+\frac{1}{2}})
    \cdot \frac{\Delta t \cdot K_i}{\Delta x} 
```

where:

- `p` is the acoustic pressure;
- `c` is the speed of sound in the medium;
- `x` is the propagation coordinate;
- `K` is the bulk modulus;
- `t` is time.

For heterogeneous domains, the acoustic properties vary spatially. The script therefore assigns each point of the domain to a material layer with its own speed of sound, density, acoustic impedance, and attenuation-related parameters, depending on the level of detail used in the simulation.

The user can place several proves across the domains so that it can track in the time domain the pressure variation.

The model is based on the methodological framework described in:

> Fernandes, N.A.T.C., Arieira, A., Hinckel, B., Silva, F., Leal, A., Carvalho, Ó. (2026). A Python FDTD Method Algorithm for 1D Planar Acoustic Wave Propagation: Simulating High-Frequency Ultrasound in the Brain and Beyond. In: Dimitrovová, Z., Biswas, P., Silva, T.A.N. (eds) *Proceedings of ICOVP and WMVC 2025*. WMVC 2025. Mechanisms and Machine Science, vol 197. Springer, Cham. https://doi.org/10.1007/978-3-032-13225-3_39

## Material properties

The material parameters used in the script are separated into two groups.

### PZT transducer parameters

The values used to represent the PZT transducer were taken from literature sources. These parameters may include, depending on the implementation:

- density;
- longitudinal wave speed;
- acoustic impedance;
- excitation frequency;
- transducer thickness or representative propagation length;
- imposed pressure or velocity amplitude.

These values should be interpreted as representative literature-based parameters rather than manufacturer-calibrated values for a specific commercial transducer. For device-specific studies, users should replace them with datasheet values or experimentally measured properties.

### Biological tissue parameters

The values used for biological tissues were taken from:

> Fernandes, N.A.T.C., Arieira, A., Leal, A., Silva, F., Carvalho, Ó. (2025). Experimental Validation of Time-Explicit Ultrasound Attenuation Modelling in Biological Tissues: A Comparative Study of Sound Diffusivity and Viscous Dissipation. *Bioengineering*, 12(9), 946. https://doi.org/10.3390/bioengineering12090946

The biological tissue parameters may include:

- density;
- speed of sound;
- attenuation coefficient.

The tissue values are intended to support ultrasound propagation studies in biological media and tissue-mimicking modelling contexts.

## How the script works

The script follows a typical FDTD workflow.

First, the spatial domain is defined. The user specifies the total propagation length and spatial resolution. The domain is then divided into discrete grid points.

Second, the material distribution is created. Each layer is assigned a thickness and acoustic properties. The script maps these layers onto the spatial grid, producing spatially varying arrays for the relevant acoustic parameters.

Third, the time step is calculated. The stability of an explicit FDTD simulation depends on the Courant-Friedrichs-Lewy condition. In one-dimensional acoustic simulations, this is commonly controlled using:

```math
\Delta t \leq \frac{\Delta x}{c_{\max}},
```

where:

- `Δt` is the time step;
- `Δx` is the spatial step;
- `c_max` is the maximum speed of sound in the simulated domain.

Fourth, the source signal is generated. The excitation may be implemented as a sinusoidal burst, continuous harmonic signal, Gaussian pulse, Ricker wavelet, or user-defined pressure/velocity input.

Fifth, the pressure field is updated over time. At each time step, the pressure at each grid point is calculated from neighbouring grid points using finite-difference approximations of the acoustic wave equation.

Sixth, boundary conditions are applied. Depending on the script version, the boundaries may be simplified as reflective, absorbing, or damped regions. Absorbing conditions are useful to reduce artificial reflections from the edges of the numerical domain.

Finally, the script stores and visualizes the results. Typical outputs include:

- pressure as a function of space and time;
- pressure measured at selected receiver positions;
- transmitted signal after each tissue layer;
- time-of-flight;
- amplitude reduction;
- frequency-domain response using FFT;
- comparison between input and received signals.

## Typical outputs

The script can be used to obtain:

- spatial pressure maps;
- temporal pressure signals at virtual sensor positions;
- propagation delay through layered tissues;
- pressure transmission across interfaces;
- attenuation trends as a function of tissue thickness;
- comparison between different material stacks;
- preliminary estimation of acoustic energy reaching a target region.

These outputs are useful for understanding how tissue properties and geometry influence ultrasound propagation before moving to more complex two-dimensional, three-dimensional, or patient-specific models.

## Installation

Clone the repository:

```bash
git clone https://github.com/your-username/your-repository-name.git
cd your-repository-name
```

Install the required Python packages:

```bash
pip install numpy scipy matplotlib pandas
```

Depending on the script version, additional packages may be required.

## Usage

Run the main script using:

```bash
python main.py
```

or, if the code is provided as a notebook:

```bash
jupyter notebook
```

Then open the notebook and execute the cells sequentially.

Before running the simulation, check the section of the code where the following parameters are defined:

- spatial resolution;
- total simulation time;
- excitation frequency;
- source amplitude;
- material layer thicknesses;
- material acoustic properties;
- receiver positions;
- plotting options.

## Example workflow

A typical simulation workflow is:

1. Define the excitation signal representing the PZT transducer.
2. Define the layered propagation path.
3. Assign acoustic properties to each layer.
4. Run the FDTD simulation.
5. Extract pressure at the target tissue or receiver position.
6. Compare the transmitted signal with the input signal.
7. Repeat the simulation for alternative tissue stacks, frequencies, or transducer amplitudes.

## Notes on interpretation

This model is a one-dimensional planar approximation. It assumes that the dominant propagation direction is normal to the tissue layers. Therefore, it does not capture:

- beam spreading;
- diffraction;
- focusing;
- complex anatomical curvature;
- shear wave conversion in solids;
- full three-dimensional scattering;
- patient-specific geometry;
- nonlinear thermal or bioheat effects, unless explicitly added.

For this reason, the model should be interpreted as a rapid computational screening tool rather than a replacement for full three-dimensional acoustic-structure interaction simulations or experimental validation.

## Recommended citation

If you use this code, method, or adapted versions of this script, please cite:
```bibtex
@software{fernandes_1D_FDTD_Multilayer,
  author       = {Fernandes, Nuno A. T. C.},
  title        = {1D FDTD Planar Acoustic Wave Propagation},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20445080},
  url          = {https://doi.org/10.5281/zenodo.20445080}
}

@incollection{Fernandes2026FDTD,
  author    = {Fernandes, Nuno A. T. C. and Arieira, A. and Hinckel, B. and Silva, F. and Leal, A. and Carvalho, {\'O}.},
  title     = {A Python FDTD Method Algorithm for 1D Planar Acoustic Wave Propagation: Simulating High-Frequency Ultrasound in the Brain and Beyond},
  booktitle = {Proceedings of ICOVP and WMVC 2025},
  editor    = {Dimitrovov{\'a}, Z. and Biswas, P. and Silva, T. A. N.},
  series    = {Mechanisms and Machine Science},
  volume    = {197},
  publisher = {Springer},
  address   = {Cham},
  year      = {2026},
  doi       = {10.1007/978-3-032-13225-3_39}
}
```


```bibtex
@incollection{Fernandes2026FDTD,
  author    = {Fernandes, Nuno A. T. C. and Arieira, A. and Hinckel, B. and Silva, F. and Leal, A. and Carvalho, {\'O}.},
  title     = {A Python FDTD Method Algorithm for 1D Planar Acoustic Wave Propagation: Simulating High-Frequency Ultrasound in the Brain and Beyond},
  booktitle = {Proceedings of ICOVP and WMVC 2025},
  editor    = {Dimitrovov{\'a}, Z. and Biswas, P. and Silva, T. A. N.},
  series    = {Mechanisms and Machine Science},
  volume    = {197},
  publisher = {Springer},
  address   = {Cham},
  year      = {2026},
  doi       = {10.1007/978-3-032-13225-3_39}
}
```

For biological tissue acoustic properties, please also cite:

```bibtex
@article{Fernandes2025Bioengineering,
  author  = {Fernandes, Nuno A. T. C. and Arieira, A. and Leal, A. and Silva, F. and Carvalho, {\'O}.},
  title   = {Experimental Validation of Time-Explicit Ultrasound Attenuation Modelling in Biological Tissues: A Comparative Study of Sound Diffusivity and Viscous Dissipation},
  journal = {Bioengineering},
  volume  = {12},
  number  = {9},
  pages   = {946},
  year    = {2025},
  doi     = {10.3390/bioengineering12090946}
}
```

## Suggested repository structure

```text
.
├── README.md
├── main.py
├── requirements.txt
├── data/
│   └── material_properties.csv
├── figures/
│   └── example_output.png
└── examples/
    └── example_layered_tissue_simulation.py
```

## Suggested `requirements.txt`

```text
numpy
scipy
matplotlib
pandas
```

## Disclaimer

This software is provided for research and educational purposes. It is not a certified medical-device simulation tool and should not be used for clinical decision-making without appropriate verification, validation, and regulatory assessment.
