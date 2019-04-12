from Angle import AngleClass

# This is a file that illustrates the use of the Angle class
# Lattice parameters and angles of the unit cell must be provided.

# first: the lattice parameters of the calculated structure
ath = 10.456788
btha = 10.456788
cth = 11.769562
#TODO: what does this parameter do?
#Please describe it
n = 0.0
alphath = 90.000000
betath = 90.000000
gammath = 90.000000

# second: the lattice parameters of the experimental structure
aex = 10.56720
bexa = 10.56720
cex = 11.9841
alphaex = 90.00
betaex = 90.000000
gammaex = 90.00

# TODO: Description is needed!
# What does this paper do?
ascal = 1.0

# Then, files including Ustar parameters have to be provided
# The matrix elements have to be given in the following order: E11 E22 E33 E23 E13 E12.

filename_Ustar_experiment = 'inpex.txt'
filename_Ustar_theory = 'inpth.txt'

angle = AngleClass(ath=ath, btha=btha, cth=cth, alphath=alphath, betath=betath, gammath=gammath, aex=aex, bexa=bexa,
                   cex=cex, alphaex=alphaex, betaex=betaex, gammaex=gammaex, n=n, ascal=ascal,
                   filename_Ustar_experiment=filename_Ustar_experiment, filename_Ustar_theory=filename_Ustar_theory)

# This will print an output file
# TODO: please describe the file in the comment part of the function to generate it!
angle.print_outputfile("out.txt")

# This will print files including main axis components
angle.print_main_axis_components_file("outexp.txt", option="experiment")
angle.print_main_axis_components_file("outth.txt", option="theory")
