# Angle
A new tool for validating theoretically derived anisotropic displacement parameters with experiment.

This python program is used to determine the angles and therefore the agreement between the eigenvectors of theoretical and experimental thermal ellipsoids. The ADPs (anisotropic displacement parameters) are an important quantity during the refinement of crystal structures. 
For this you will need DFT-programs, such as [VASP](https://www.vasp.at/) and a phonon-code, such as [Phonopy](https://github.com/atztogo/phonopy).

# What to cite
Please cite this repository: Damian Mroz, Janine George,& Richard Dronskowski. (2019, April 08). Angle (Version 1.0).
If used, please also cite [VASP](https://www.vasp.at/) and [Phonopy](https://github.com/atztogo/phonopy).
In addition to this, you should cite [MolecularToolbox](https://github.com/JaGeo/MolecularToolbox) if used. Please follow the citation guidelines on this repository.

# Installation
You will need to export the python path of the folder including the "Angle.py" file or copy it to the folder where you run the script. Furthermore, you will need `numpy` and `math`.

# How to
You have to perform DFT calculations for example with [VASP](https://www.vasp.at/) first, followed by phonon evaluation (e.g. [Phonopy](https://github.com/atztogo/phonopy)). Furthermore, you need to transform the thermal displacement matrices ([Phonopy](https://github.com/atztogo/phonopy) output) to the "Ustar" format. You can use the [MolecularToolbox](https://github.com/JaGeo/MolecularToolbox) to obtain the Ustar values directly from [Phonopy](https://github.com/atztogo/phonopy).
Then you just need to prepare the input files "inpex.txt" and "inpth.txt" as shown in the `Example` folder. These files include the matrix components of Ustar which can be derived with the help of the [MolecularToolbox](https://github.com/JaGeo/MolecularToolbox). Please make sure that the Ustar parameters are sorted in the correct way: E11 E22 E33 E23 E13 E12.

# Result
You will get three output files: Two of them will give you the main axis components of the ADPs. The other one will give you a quantity related to the length of the semi principal axes and angles between the theoretical and experimental main axes.

# Information about the Authors
* D. Mroz (RWTH Aachen University)
* Refactoring and changes to the code by J. George (Université Catholique de Louvain)
* PI during the development of the code: R. Dronskowski, RWTH Aachen University


