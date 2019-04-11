import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv
from numpy import dot, arccos, clip
from numpy.linalg import norm
import math


class AngleClass:

    def __init__(self, aex, bexa, cex, alphaex, betaex, gammaex, ath, btha, cth, alphath, betath, gammath, n,
                 filename_Ustar_experiment, filename_Ustar_theory, ascal=1.0):
        """

        :param aex: a parameter experimental strucutre
        :param bexa: b parameter experimental structure
        :param cex: c parameter experimental structure
        :param alphaex: alpha parameter experimental structure
        :param betaex: beta parameter experimental strucutre
        :param gammaex: gamma parameter experimental structure
        :param ath: a parameter computed structure
        :param btha: b parameter computed structure
        :param cth: c parameter computed structure
        :param alphath: alpha parameter computed structure
        :param betath: beta parameter computed structure
        :param gammath: gamma paramter computed structure
        :param n: TODO: fill this in
        :param filename_Ustar_experiment: filename of Ustar parameters of experimental ADPs
        :param filename_Ustar_theory: filename of Ustar parameters of theoretical ADPs
        :param ascal: TODO: fill this in
        """
        number_atoms_theory = self._get_number_atoms(filename_Ustar_theory)
        number_atoms_experiment = self._get_number_atoms(filename_Ustar_experiment)

        if number_atoms_experiment != number_atoms_theory:
            raise ValueError("The number of ADPs in your input files is different")

        # loop over all ADPs
        self.result_dicts = []
        for iatom in range(0, number_atoms_theory):
            # get names of atoms from input files
            atomnameth = self._get_atom_name(filename=filename_Ustar_theory, atomnumber=iatom)
            atomnameex = self._get_atom_name(filename=filename_Ustar_experiment, atomnumber=iatom)

            # read matrix from file
            Aex = self._get_Umatrix_from_file(filename=filename_Ustar_experiment, atomnumber=iatom)
            Ath = self._get_Umatrix_from_file(filename=filename_Ustar_theory, atomnumber=iatom, ascal=ascal)

            # get Ucart:
            Aex22 = self._calculate_Ucart(ath=aex, btha=bexa, cth=cex, alphath=alphaex, betath=betaex, gammath=gammaex,
                                          n=n,
                                          Ath=Aex)
            Ath22 = self._calculate_Ucart(ath=ath, btha=btha, cth=cth, alphath=alphath, betath=betath, gammath=gammath,
                                          n=n,
                                          Ath=Ath)

            # get sorted eigenvalues from Ucart matrix
            wexpv = LA.eigvals(Aex22)  # Eigenvalues are calculated.
            wexps = np.sort(wexpv)  # Eigenvalues are sorted in an ascending fashion.
            wthpv = LA.eigvals(Ath22)
            wthps = np.sort(wthpv)

            # maximum value of the experimental eigenvalues is determined
            Aexmin = np.amin([[wexpv[0], wexpv[1], wexpv[2]]])
            # minimum value of the experimental eigenvalues is determined
            Aexmax = np.amax([[wexpv[0], wexpv[1], wexpv[2]]])  # Minimum value is determined.

            # get eigenvalues and eigenvectors of inverted matrix
            # x^T U_cart^{-1] x=1
            wex, vex = self._get_eigenvalues_eigenvetors_inverted_matrix(Ath22=Aex22)
            wth, vth = self._get_eigenvalues_eigenvetors_inverted_matrix(Ath22=Ath22)

            # angles between all main axes of caculated and theoretical ADP are calculated
            angle1, angle2, angle3, angle4, angle5, angle6, angle7, angle8, angle9 = self._get_angles(vth=vth, vex=vex)

            # calculates theoretical and experimental volume of the ellipsoids
            Volth = self._get_volume_ellipsoid(wth=wth)
            Volex = self._get_volume_ellipsoid(wth=wex)

            # calculate ratio of the volumes
            Ratio_Volumes = Volth / Volex

            # transforms the computed eigenvectors from the cartesian to the crystallographic coordination system
            Vectorth = self._get_transformation_to_crystallographic_CS(ath=ath, btha=btha, cth=cth, alphath=alphath,
                                                                       betath=betath, gammath=gammath)
            Vectorex = self._get_transformation_to_crystallographic_CS(ath=aex, btha=bexa, cth=cex, alphath=alphaex,
                                                                       betath=betaex, gammath=gammaex)

            # get vectors in crystallographic coordinate system
            new_cexp1, new_cexp2, new_cexp3 = self.get_vectors_in_crystal_cs(Vectorth=Vectorex, vth=vex)
            new_cth1, new_cth2, new_cth3 = self.get_vectors_in_crystal_cs(Vectorth=Vectorth, vth=vth)

            # save result for printing or other use
            self.result_dicts.append(
                {"Aexmax": Aexmax, "Aexmin": Aexmin, "wth": wth, "wex": wex, "angle1": angle1, "angle2": angle2,
                 "angle3": angle3, "angle4": angle4, "angle5": angle5, "angle6": angle6, "angle7": angle7,
                 "angle8": angle8, "angle9": angle9, "Ratio": Ratio_Volumes, "new_cexp1": new_cexp1,
                 "new_cexp2": new_cexp2,
                 "new_cexp3": new_cexp3, "new_cth1": new_cth1, "new_cth2": new_cth2, "new_cth3": new_cth3,
                 "atomname": atomnameex, "atomname_theory": atomnameth, "wexps": wexps, "wthps": wthps})

    def print_main_axis_components_file(self, filename, option='theory'):
        """
        prints a file with the main axis components
        :param filename: name of the outputfile
        :param option: "theory" and "experiment" are available
        """
        with open(filename, 'w') as f:
            for dictionary in self.result_dicts:
                if option == "experiment":
                    print(dictionary["atomname"], dictionary["wexps"][0], dictionary["wexps"][1],
                          dictionary["wexps"][2], file=f)
                elif option == "theory":
                    print(dictionary["atomname_theory"], dictionary["wthps"][0], dictionary["wthps"][1],
                          dictionary["wthps"][2], file=f)
                else:
                    raise ValueError("wrong option")

    def print_outputfile(self, filename="out.txt"):
        """
        TODO: please describe the output file
        :param filename: filename of the output file
        """
        with open(filename, 'w') as f:
            for dictionary in self.result_dicts:
                if dictionary["Aexmax"] / dictionary["Aexmin"] > 1.8:
                    print("Anisotropic", dictionary["atomname"], dictionary["Aexmax"] / dictionary["Aexmin"],
                          "Experiment1", dictionary["new_cexp1"], dictionary["wex"][0], "Experiment2",
                          dictionary["new_cexp2"], dictionary["wex"][1], "Experiment3", dictionary["new_cexp3"],
                          dictionary["wex"][2], "Theory1", dictionary["new_cth1"], dictionary["wth"][0], "Theory2",
                          dictionary["new_cth2"], dictionary["wth"][1], "Theory3", dictionary["new_cth3"],
                          dictionary["wth"][2], file=f)
                print(dictionary["atomname"], "Ratio", dictionary["Ratio"], "Angle", "(11)",
                    dictionary["angle1"] * 57.296, "(22)", dictionary["angle5"] * 57.296, "(33)",
                    dictionary["angle9"] * 57.296, "(12)", dictionary["angle2"] * 57.296, "(21)",
                    dictionary["angle4"] * 57.296, "(13)", dictionary["angle3"] * 57.296, "(32)",
                    dictionary["angle8"] * 57.296, file=f)
                print("(23)", dictionary["angle6"] * 57.296, "(31)", dictionary["angle7"] * 57.296,
                    file=f)  # Angles are given in degrees; format is experiment theory (12 for example means experiment 1, theory 2). Anisotropic ellipsoids are highlighted. The vectors are transformed to the crystallographic coordinate system.

    def _get_number_atoms(self, filename):
        return sum(1 for line in open(filename) if line != '\n')

    def _get_Umatrix_from_file(self, filename, atomnumber, ascal=1.0):
        with open(filename, 'r') as f:
            for iline, line in enumerate(f):
                if iline == atomnumber:
                    splitline = line.split()
                    newarray = splitline[1:]
                    matrix = np.array(
                        [[ascal * float(newarray[0]), ascal * float(newarray[5]), ascal * float(newarray[4])],
                         [ascal * float(newarray[5]), ascal * float(newarray[1]), ascal * float(newarray[3])],
                         [ascal * float(newarray[4]), ascal * float(newarray[3]), ascal * float(newarray[2])]])
        return matrix

    def _get_atom_name(self, filename, atomnumber):
        with open(filename, 'r') as f:
            for iline, line in enumerate(f):
                if iline == atomnumber:
                    splitline = line.split()
                    atom = splitline[0]
        return atom

    def _get_transformation_matrix_to_Cartesian_CS(self, ath, btha, cth, alphath, betath, gammath, n):

        B2th = np.array([[ath, btha * (np.cos(gammath * np.pi / 180.0)), cth * np.cos(betath * np.pi / 180.0)],
                         [n, btha * np.sin(gammath * np.pi / 180.0), cth * ((
                                 np.cos(alphath * np.pi / 180.0) - (np.cos(betath * np.pi / 180.0)) * np.cos(
                             gammath * np.pi / 180.0))) / np.sin(gammath * np.pi / 180.0)], [n, n, (
                    ath * btha * cth * np.sqrt(
                1 - np.cos(alphath * np.pi / 180.0) * np.cos(alphath * np.pi / 180.0) - np.cos(
                    betath * np.pi / 180.0) * np.cos(betath * np.pi / 180.0) - np.cos(gammath * np.pi / 180.0) * np.cos(
                    gammath * np.pi / 180.0) + 2 * np.cos(alphath * np.pi / 180.0) * np.cos(
                    betath * np.pi / 180.0) * np.cos(gammath * np.pi / 180.0))) / (ath * btha * np.sin(
                gammath * np.pi / 180.0))]])  # Conversion to Cartesian coordinates.

        return B2th

    def _calculate_Ucart(self, ath, btha, cth, alphath, betath, gammath, n, Ath):
        B2th = self._get_transformation_matrix_to_Cartesian_CS(ath=ath, btha=btha, cth=cth, alphath=alphath,
                                                               betath=betath, gammath=gammath, n=n)

        Mth = np.transpose(B2th)  # The matrix is transposed.
        Ath2 = np.matmul(B2th, Ath)
        Ath22 = np.matmul(Ath2, Mth)

        return Ath22

    def _get_eigenvalues_eigenvetors_inverted_matrix(self, Ath22):
        Ath = inv(Ath22)
        wth, vth = LA.eig(Ath)
        return wth, vth

    def _reformat_eigenvectors(self, vth):
        uth1 = np.array([vth[0, 0], vth[1, 0], vth[2, 0]])
        uth2 = np.array([vth[0, 1], vth[1, 1], vth[2, 1]])
        uth3 = np.array([vth[0, 2], vth[1, 2], vth[2, 2]])
        return uth1, uth2, uth3

    def _get_angles(self, vth, vex):
        uex1, uex2, uex3 = self._reformat_eigenvectors(vex)
        uth1, uth2, uth3 = self._reformat_eigenvectors(vth)

        c1 = dot(uex1, uth1) / norm(uex1) / norm(uth1)
        c2 = dot(uex1, uth2) / norm(uex1) / norm(uth2)
        c3 = dot(uex1, uth3) / norm(uex1) / norm(uth3)
        c4 = dot(uex2, uth1) / norm(uex2) / norm(uth1)
        c5 = dot(uex2, uth2) / norm(uex2) / norm(uth2)
        c6 = dot(uex2, uth3) / norm(uex2) / norm(uth3)
        c7 = dot(uex3, uth1) / norm(uex3) / norm(uth1)
        c8 = dot(uex3, uth2) / norm(uex3) / norm(uth2)
        c9 = dot(uex3, uth3) / norm(uex3) / norm(uth3)
        # Angles are calculated with the arccos. The numbers denote all possible permutations.
        angle1 = arccos(clip(c1, -1, 1))
        angle2 = arccos(clip(c2, -1, 1))
        angle3 = arccos(clip(c3, -1, 1))
        angle4 = arccos(clip(c4, -1, 1))
        angle5 = arccos(clip(c5, -1, 1))
        angle6 = arccos(clip(c6, -1, 1))
        angle7 = arccos(clip(c7, -1, 1))
        angle8 = arccos(clip(c8, -1, 1))
        angle9 = arccos(clip(c9, -1, 1))

        return angle1, angle2, angle3, angle4, angle5, angle6, angle7, angle8, angle9

    def _get_volume_ellipsoid(self, wth):
        # Ellipsoid volume is calculated.
        Volth = (1 / (math.sqrt(abs(wth[0])))) * (1 / (math.sqrt(abs(wth[1])))) * (
                1 / (math.sqrt(abs(wth[2]))))
        return Volth

    def _get_transformation_to_crystallographic_CS(self, ath, btha, cth, alphath, betath, gammath):
        Vectorth = np.array([[1 / ath, -(np.cos(gammath * np.pi / 180.0)) / (ath * np.sin(gammath * np.pi / 180.0)), (
                btha * cth * (np.cos(alphath * np.pi / 180.0) * np.cos(gammath * np.pi / 180.0) - np.cos(
            betath * np.pi / 180.0))) / (np.sin(gammath * np.pi / 180.0) * ath * btha * cth * np.sqrt(
            1 - np.cos(alphath * np.pi / 180.0) * np.cos(alphath * np.pi / 180.0) - np.cos(
                betath * np.pi / 180.0) * np.cos(
                betath * np.pi / 180.0) - np.cos(gammath * np.pi / 180.0) * np.cos(
                gammath * np.pi / 180.0) + 2 * np.cos(
                alphath * np.pi / 180.0) * np.cos(betath * np.pi / 180.0) * np.cos(gammath * np.pi / 180.0)))],
                             [0, 1 / (btha * np.sin(gammath * np.pi / 180.0)), ath * cth * (
                                     np.cos(betath * np.pi / 180.0) * np.cos(gammath * np.pi / 180.0) - np.cos(
                                 alphath * np.pi / 180.0)) / (
                                      np.sin(gammath * np.pi / 180.0) * ath * btha * cth * np.sqrt(
                                  1 - np.cos(alphath * np.pi / 180.0) * np.cos(alphath * np.pi / 180.0) - np.cos(
                                      betath * np.pi / 180.0) * np.cos(betath * np.pi / 180.0) - np.cos(
                                      gammath * np.pi / 180.0) * np.cos(gammath * np.pi / 180.0) + 2 * np.cos(
                                      alphath * np.pi / 180.0) * np.cos(betath * np.pi / 180.0) * np.cos(
                                      gammath * np.pi / 180.0)))], [0, 0,
                                                                    (ath * btha * np.sin(gammath * np.pi / 180.0)) / (
                                                                            ath * btha * cth * np.sqrt(1 - np.cos(
                                                                        alphath * np.pi / 180.0) * np.cos(
                                                                        alphath * np.pi / 180.0) - np.cos(
                                                                        betath * np.pi / 180.0) * np.cos(
                                                                        betath * np.pi / 180.0) - np.cos(
                                                                        gammath * np.pi / 180.0) * np.cos(
                                                                        gammath * np.pi / 180.0) + 2 * np.cos(
                                                                        alphath * np.pi / 180.0) * np.cos(
                                                                        betath * np.pi / 180.0) * np.cos(
                                                                        gammath * np.pi / 180.0)))]])

        return Vectorth

    def get_vectors_in_crystal_cs(self, Vectorth, vth):
        uth1, uth2, uth3 = self._reformat_eigenvectors(vth)
        new_vec1 = np.matmul(Vectorth, uth1)
        new_vec2 = np.matmul(Vectorth, uth2)
        new_vec3 = np.matmul(Vectorth, uth3)
        return new_vec1, new_vec2, new_vec3
