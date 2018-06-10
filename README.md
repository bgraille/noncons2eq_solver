
These python scripts can be used to retrieve the numerical results of the paper untitled 

_Numerical treatment of the nonconservative product in a multiscale fluid model for plasmas in thermal nonequilibrium: application to solar physics_

authors: Q. Wargnier, S. Faure, B. Graille, T. Magin, M. Massot

*How use it*
> * clone the repository or download it
> * `python  main.py`
> * choose a test case corresponding to a json file

The json files are used to store the parameters of the simulation with the following conventions

```
{
  "version":"1.0",                   ### version of the domcument
  "prefix":"results/",               ### where the results have to be store (pictures)
  "Lx":     10.0,                    ### length of the domain in space
  "N":      2000,                    ### number of points in space
  "t0":     0.0,                     ### initial time
  "tf":     1.0,                     ### final time
  "peraff": 0.1,                     ### time laps between plots
  "gamma":  1.6666666666666667,      ### value of the adiabatic coefficient gamma = 5/3
  "diffD":  0.1,                     ### value of the electron diffusion coefficient
  "diffL":  0.001,                   ### value of the electron thermal conductivity
  "onde":   3,                       ### type of the wave
  "rhoeL":  "None",                  ### value of the electron mass (left)
  "rhohL":  "None",                  ### value of the heavy particle mass (left)
  "peL":    "None",                  ### value of the electron pressure (left)
  "pL":     1.5,                     ### value of the total pressure (left)
  "vhL":    "None",                  ### value of the hydrodynamic velocity (left)
  "rhoeR":  0.01,                    ### value of the electron mass (right)
  "rhohR":  1.0,                     ### value of the heavy particle mass (right)
  "peR":    0.1,                     ### value of the electron pressure (right)
  "pR":     1.0,                     ### value of the total pressure (right)
  "vhR":    0.2,                     ### value of the hydrodynamic velocity (right)
  "splitting":     "True",           ### use operator splitting
  "diff_substeps": "True",           ### use time substeps for the diffusive part (only with splitting)
  "scheme_conv":   "LF",             ### choice of the convective scheme (LF: Lax-Friedrichs, UW: Upwind)
  "scheme_nc":     "df",             ### choice of the diffusive scheme (df: classical finite differences, MW: our scheme)
  "CFL":    0.1,                     ### link between time and space step (the scale depends of the splitting)
  "liste_plot": ["rhoe", "pe", "Te"] ### list of the outputs
}
```