
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
  "version":"1.0",
  "prefix":"results/",
  "peraff": 0.1,
  "Lx":     10.0,
  "N":      2000,
  "t0":     0.0,
  "tf":     1.0,
  "CFL":    0.1,
  "gamma":  1.6666666666666667,
  "diffD":  0.1,
  "diffL":  0.001,
  "onde":   3,
  "rhoeL":  "None",
  "rhohL":  "None",
  "peL":    "None",
  "pL":     1.5,
  "vhL":    "None",
  "rhoeR":  0.01,
  "rhohR":  1.0,
  "peR":    0.1,
  "pR":     1.0,
  "vhR":    0.2,
  "diff_substeps": "True",
  "splitting":     "True",
  "scheme_conv":   "LF",
  "scheme_nc":     "df",
  "liste_plot": ["rhoe", "pe", "Te"]
}
```