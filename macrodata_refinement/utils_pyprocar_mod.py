
class pyprocar:
    '''
    The whole script for reading and processing the PROCAR file is appart, in this class.
    Great part of it is based in the pyprocar code from 5/12/2014
    '''
    def __init__(self):
        # self.fileStr = './PROCAR'
        pass

    def _read_kpoints(self, permissive):
        import re
        """
        Reads the k-point headers. A typical k-point line is:
        k-point    1 :    0.00000000 0.00000000 0.00000000  weight = 0.00003704\n

        fills self.kpoint[kpointsCount][3]
        
        The weights are discarded (are they useful?)
        """
        # finding all the K-points headers 
        self.kpoints = 0
        self.kpoints = re.findall(r"k-point\s+\d+\s*:\s+([-.\d\s]+)", self.fileStr)
        # self.log.debug(str(len(self.kpoints))+" K-point headers found")
        # self.log.debug("The first match found is: " + str(self.kpoints[0]))
        # trying to build an array
        self.kpoints = [x.split() for x in self.kpoints]
        
        try:
            self.kpoints = np.array(self.kpoints, dtype=float)
        except ValueError:
            # self.log.error("Ill-formatted data:")
            print("\n".join([str(x) for x in self.kpoints]))
            if permissive is True:
                # Discarding the kpoints list, however I need to set
                # self.ispin beforehand.
                if len(self.kpoints) == self.kpointsCount:
                    self.ispin = 1
                elif len(self.kpoints) == 2 * self.kpointsCount:
                    self.ispin = 2
                else:
                    raise ValueError("Kpoints do not match with ispin=1 or 2.")
                self.kpoints = None
                self.log.warning("K-points list is useless, setting it to `None`")
                return
            else:
                raise ValueError("Badly formated Kpoints headers, try `--permissive`")
        # # if successful, go on

        # trying to identify an non-polarized or non-collinear case, a
        # polarized case or a defective file

        if len(self.kpoints) != self.kpointsCount:
            # if they do not match, may means two things a spin polarized
            # case or a bad file, lets check
            # self.log.debug(
            #     "Number of kpoints do not match, looking for a " "spin-polarized case"
            # )
            # lets start testing if it is spin polarized, if so, there
            # should be 2 identical blocks of kpoints.
            up, down = np.vsplit(self.kpoints, 2)
            if (up == down).all():
                # self.log.info("Spin-polarized calculation found")
                self.ispin = 2
                # just keeping one set of kpoints (the other will be
                # discarded)
                self.kpoints = up
            else:
                # self.log.error("Number of K-points do not match! check them.")
                raise RuntimeError("Bad Kpoints list.")
        # if ISPIN != 2 setting ISPIN=1 (later for the non-collinear case 1->4)
        # It is unknown until parsing the projected data
        else:
            self.ispin = 1
        #checking again, for compatibility,
        if len(self.kpoints) != self.kpointsCount:
            raise RuntimeError(
                "Kpoints number do not match with metadata (header of PROCAR)"
            )

        # self.log.debug(str(self.kpoints))
        # self.log.info("The kpoints shape is " + str(self.kpoints.shape))

        # if self.recLattice is not None:
        #     # self.log.info("Changing to cartesians coordinates")
        #     self.kpoints = np.dot(self.kpoints, self.recLattice)
        #     # self.log.debug("New kpoints: \n" + str(self.kpoints))
        return

    def _read_bands(self):
        """Reads the bands header. A typical bands is:
        band   1 # energy   -7.11986315 # occ.  1.00000000

        fills self.bands[kpointsCount][bandsCount]

        The occupation numbers are discarded (are they useful?)"""
        # self.log.debug("readBands")
        # if not self.fileStr:
        #     log.warning("You should invoke `procar.read()` instead. Returning")
        #     return

        # finding all bands
        self.bands = re.findall(
            r"band\s*(\d+)\s*#\s*energy\s*([-.\d\s]+)", self.fileStr
        )
        # self.log.debug(
        #     str(len(self.bands))
        #     + " bands headers found, bands*Kpoints = "
        #     + str(self.bandsCount * self.kpointsCount)
        # )
        # self.log.debug("The first match found is: " + str(self.bands[0]))

        # checking if the number of bands match

        if len(self.bands) != self.bandsCount * self.kpointsCount * self.ispin:
            # self.log.error("Number of bands headers do not match")
            raise RuntimeError("Number of bands don't match")

        # casting to array to manipulate the bands
        self.bands = np.array(self.bands, dtype=float)
        # self.log.debug(str(self.bands))

        # Now I will deal with the spin polarized case. The goal is join
        # them like for a non-magnetic case
        if self.ispin == 2:
            # up and down are along the first axis
            up, down = np.vsplit(self.bands, 2)
            # self.log.debug("up   , " + str(up.shape))
            # self.log.debug("down , " + str(down.shape))

            # reshapping (the 2  means both band index and energy)
            up.shape = (self.kpointsCount, self.bandsCount, 2)
            down.shape = (self.kpointsCount, self.bandsCount, 2)

            # setting the correct number of bands (up+down)
            self.bandsCount *= 2
            # self.log.debug("New number of bands : " + str(self.bandsCount))

            # and joining along the second axis (axis=1), ie: bands-like
            self.bands = np.concatenate((up, down), axis=1)

        # otherwise just reshaping is needed
        else:
            self.bands.shape = (self.kpointsCount, self.bandsCount, 2)

        # Making a test if the broadcast is rigth, otherwise just print
        test = [x.max() - x.min() for x in self.bands[:, :, 0].transpose()]
        # if np.array(test).any():
            # self.log.warning(
            #     "The indexes of bands do not match. CHECK IT. "
            #     "Likely the data was wrongly broadcasted"
            # )
            # self.log.warning(str(self.bands[:, :, 0]))
        # Now safely removing the band index
        self.bands = self.bands[:, :, 1]
        # self.log.info("The bands shape is " + str(self.bands.shape))
        return

    def _read_orbitals(self):
        self.orbitalName = [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "x2-y2",
            "fy3x2",
            "fxyz",
            "fyz2",
            "fz3",
            "fxz2",
            "fzx2",
            "fx3",
            "tot",
        ]
        self.orbitalName_old = [
            "s",
            "py",
            "pz",
            "px",
            "dxy",
            "dyz",
            "dz2",
            "dxz",
            "dx2",
            "tot",
        ]
        self.orbitalName_short = ["s", "p", "d", "f", "tot"]
        """Reads all the spd-projected data. A typical/expected block is:
        ion      s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot
        1  0.079  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.079
        2  0.152  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.152
        3  0.079  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.079
        4  0.188  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.188
        5  0.188  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.188
        tot  0.686  0.000  0.002  0.000  0.000  0.000  0.000  0.000  0.000  0.688
        (x2 for spin-polarized -akwardkly formatted-, x4 non-collinear -nicely
        formatted-).

        The data is stored in an array self.spd[kpoint][band][ispin][atom][orbital]

        Undefined behavior in case of phase factors (LORBIT = 12).
        """
        # finding all orbital headers
        self.spd = re.findall(r"ion(.+)", self.fileStr)

        # testing if the orbital names are known (the standard ones)
        FoundOrbs = self.spd[0].split()
        size = len(FoundOrbs)
        # only the first 'size' orbital
        StdOrbs = self.orbitalName[: size - 1] + self.orbitalName[-1:]
        StdOrbs_short = self.orbitalName_short[: size - 1] + self.orbitalName_short[-1:]
        StdOrbs_old = self.orbitalName_old[: size - 1] + self.orbitalName_old[-1:]

        self.orbitalCount = size
        self.orbitalNames = self.spd[0].split()
        # Now reading the bulk of data
        # The case of just one atom is handled differently since the VASP
        # output is a little different
        if self.ionsCount == 1:
            self.spd = re.findall(r"^(\s*1\s+.+)$", self.fileStr, re.MULTILINE)
        else:
            self.spd = re.findall(r"([-.\d\se]+tot.+)\n", self.fileStr)
        # free the memory (could be a lot)
        self.fileStr = None
        expected = self.bandsCount*self.kpointsCount
        if expected == len(self.spd):
            pass
        elif expected * 4 == len(self.spd):
            self.ispin = 4
        else:
            raise RuntimeError("Shit happens")

        # checking for consistency
        for line in self.spd:
            if len(line.split()) != (self.ionsCount) * (self.orbitalCount + 1):
                print(line)
                raise RuntimeError("Flats happens")

        # replacing the "tot" string by a number, to allows a conversion
        # to numpy
        self.spd = [x.replace("tot", "0") for x in self.spd]
        self.spd = [x.split() for x in self.spd]
        self.spd = np.array(self.spd, dtype=float)
        # self.log.debug("The spd (old) array shape is:" + str(self.spd.shape))

        # handling collinear polarized case
        if self.ispin == 2:
            # self.log.debug("Handling spin-polarized collinear case...")
            # splitting both spin components, now they are along k-points
            # axis (1st axis) but, then should be concatenated along the
            # bands.
            up, down = np.vsplit(self.spd, 2)
            # ispin = 1 for a while, we will made the distinction
            up.shape = (self.kpointsCount, self.bandsCount/2, 1,
                        self.ionsCount, self.orbitalCount+1)
            down.shape = (self.kpointsCount, self.bandsCount/2, 1,
                        self.ionsCount, self.orbitalCount+1)

            density = np.concatenate((up, down), axis=1)
            magnet = np.concatenate((up, -down), axis=1)
            # concatenated along 'ispin axis'
            self.spd = np.concatenate((density, magnet), axis=2)

        # otherwise, just a reshaping suffices
        else:
            self.spd.shape=(self.kpointsCount, self.bandsCount, self.ispin,
                      self.ionsCount, self.orbitalCount+1)

        # self.log.info("spd array ready. Its shape is:" + str(self.spd.shape))
        return

    def _read_spin(self):
        pass

    def _read_file(self,
        path_read,
        permissive=False):

        file = open(path_read+'/PROCAR')
        # Line 1: PROCAR lm decomposed
        file.readline()  # throwaway
        # Line 2: # of k-points:  816   # of bands:  52   # of ions:   8
        metaLine = file.readline() # metadata
        # self.log.debug("The metadata line is: "+ metaLine)
        re.findall(r"#[^:]+:([^#]+)", metaLine)
        self.kpointsCount, self.bandsCount, self.ionsCount = map(
            int, re.findall(r"#[^:]+:([^#]+)", metaLine))
        
        # self.log.info("kpointsCount = " + str(self.kpointsCount));
        # self.log.info("bandsCount = " + str(self.bandsCount));
        # self.log.info("ionsCount = " + str(self.ionsCount));
        if self.ionsCount == 1:
            pass
            # self.log.warning("Special case: only one atom found. The program may not work as expected")
        else:
            # self.log.debug("An extra ion will be the total value")
            self.ionsCount = self.ionsCount + 1

        # #reading all the rest of the file to be parsed below
        self.fileStr = file.read()
        self._read_kpoints(permissive)
        self._read_bands()
        self._read_orbitals()
        # self.log.debug("readfile...done")
        return

    def filter_procar(self):
        pass

    def process_procar(self):
        pass

