# Angle
A new tool for validating theoretically derived anisotropic displacement parameters with experiment.

This python program is used to determine the angles and therefore the agreement between the eigenvectors of theoretical and experimental thermal ellipsoids. The ADPs (anisotropic displacement parameters) are an important quantity during the refinement of crystal structures. 
For this you will need DFT-programs, such as [VASP](https://www.vasp.at/) and a phonon-code, such as [Phonopy](https://github.com/atztogo/phonopy).

# What to cite
Please cite this repsoitory: Damian Mroz, & Richard Dronskowski. (2019, April 08). Angle (Version 1.0). 
If used, please also cite [VASP](https://www.vasp.at/) and [Phonopy](https://github.com/atztogo/phonopy).
In addition to this, you should cite [MolecularToolbox](https://github.com/JaGeo/MolecularToolbox) if used. In this case please cite: Janine George, & Richard Dronskowski. (2017, November 30). Molecular Toolbox (Version 1.0.2). Zenodo. [http://doi.org/10.5281/zenodo.1069052](http://doi.org/10.5281/zenodo.1069052) (BibTeX).

# Installation
You can directly run the python script, but you will need `numpy`.

# How to
You have to perform DFT calculations for example with [VASP](https://www.vasp.at/) first, followed by phonon evaluation (e.g. [Phonopy](https://github.com/atztogo/phonopy). Furthermore, you need to transform the thermal displacement matrices ([Phonopy](https://github.com/atztogo/phonopy) output) to the "Ustar" format. You can use the [MolecularToolbox](https://github.com/JaGeo/MolecularToolbox) to obtain the Ustar values directly from [Phonopy](https://github.com/atztogo/phonopy). Then you just need to prepare the input filex "inpex.txt" and "inpth.txt" as shown in the Examples.

# Result
You will get a list with the eigenvectors and corresponding angles.

# Information about the Author
 *D. Mroz (RWTH Aachen University)
*PI during the development of the code: R. Dronskowski, RWTH Aachen University


