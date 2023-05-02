# GrAED
Grid Adsorption Energy Descriptors (GrAED) is an adsorption energy sampling tool based on symmetric grids to calculate the adsorption enthalpy and the Henry constants in nanoporous materials as well as statistical quantities used in ML modelling of adsorption process

## Installation

Check if you have `c++11` compiler installed (may work with other compiler but has been tested and mainly used with this compiler).

Tested for c++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0

Compilation that follows the rules set in the `Makefile`:
```
make all
```

## Example Usage

If you want to run a surface sampling simulation on the structure KAXQIL (CSD code) from CoRE MOF 2019 all-solvent removed with the Dreiding+uff forcefield at 298K with a 12A cutoff for the xenon.
```
./graed structure/KAXQIL_clean_14.cif forcefield/UFF.def 298 12.0 Xe 0.12 100 0.8
```
The inputs correspond to the structure cif file, the forcefield definition file (Raspa), the temperature, the cutoff in LJ potential, the atom, the grid spacing, the energy threshold under which the nergy is accurately calculated, the relative size of the blocking framework atoms (between 0 and 1).

You should get an output that has values close to this:
```
KAXQIL_clean_14,-44.626,2.21797,2.14037,10.0172,15.0287,73.1689,0.0300725,0.218192
```
The results are printed in a comma separated format: structure name, adsorption enthalpy (kJ/mol), standard deviation of Boltzmann weighted energies, skewness, kurtosis, standard deviation of energies, Henry coefficient (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)

This binary is supposed to be used in a high-throughput manner to add rows to a csv file.

## Acknowledgement

This code has been developped during a PhD thesis co-financed by the CEA and Orano under the supervision of François-Xavier Coudert: https://github.com/fxcoudert

This code includes the library developped in Gemmi:
https://github.com/project-gemmi/gemmi.git

## License

The MIT License (MIT)

Copyright (c) 2023 Emmanuel Ren

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
