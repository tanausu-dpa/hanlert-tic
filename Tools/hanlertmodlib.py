import struct,copy
import numpy as np

######################################################################
######################################################################
######################################################################
######################################################################

def airtovacuum(wave,limit=None):
    ''' Convert from air to vacuum, VALD3
        http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    '''

    if limit is None:
        limit = 2000.0

    s = 1e4/wave
    s2 = s*s
    ff = 1.0 + 8.336624212083e-5 + \
               2.408926869968e-2/(130.10659245522 - s2) + \
               1.599740894897e-4/(38.92568793293 - s2)
    out = copy.deepcopy(wave)

    for i in range(len(wave)):
        if out[i] > limit:
            out[i] *= ff[i]

    return out

######################################################################
######################################################################
######################################################################
######################################################################

class allen_class():
    ''' Compute specific intensity from Allen
    '''

    def __init__(self):
        ''' Create interpolation functions to be used later.
            Inputs: None
            Outputs: None
            Internal: Defines interpolation functions for CLV
                      coefficients and continuum intensity (in
                      wavelength)
        '''

        import scipy.interpolate as inter

        # Size
        N = 22

        # Constants
        self.c = 2.99792458e10  # speed of light [cm/s]

        # Lambda [micron]
        xI = np.array([0.2e0,0.22e0,0.24e0,0.26e0,0.28e0,0.3e0, \
                       0.32e0,0.34e0,0.36e0,0.37e0,0.38e0,0.39e0, \
                       0.4e0,0.41e0,0.42e0,0.43e0,0.44e0,0.45e0, \
                       0.46e0,0.48e0,0.5e0,0.55e0,0.6e0,0.65e0, \
                       0.7e0,0.75e0,0.8e0,0.9e0,1e0,1.1e0,1.2e0, \
                       1.4e0,1.6e0,1.8e0,2e0,2.5e0,3e0,4e0,5e0, \
                       6e0,8e0,10e0,12e0])

        # Intensity [10^10 erg s^-1 cm^-2 sr^-1 micron^-1] TODO -> {CHECK}
        yI = np.array([0.06e0,0.21e0,0.29e0,0.6e0,1.3e0,2.45e0, \
                       3.25e0,3.77e0,4.13e0,4.23e0,4.63e0,4.95e0, \
                       5.15e0,5.26e0,5.28e0,5.24e0,5.19e0,5.1e0, \
                       5e0,4.79e0,4.55e0,4.02e0,3.52e0,3.06e0, \
                       2.69e0,2.28e0,2.03e0,1.57e0,1.26e0,1.01e0, \
                       0.81e0,0.53e0,0.36e0,0.238e0,0.16e0,0.078e0, \
                       0.041e0,0.0142e0,0.0062e0,0.0032e0,0.00095e0, \
                       0.00035e0,0.00018e0])

        # Convert units to cgs in Hz bandwidth
        yI *= 1e10*xI*xI/(self.c*1e4)

        # Interpolation functions
        self.fI = inter.interp1d(xI,yI,kind='linear', \
                        bounds_error=False, \
                        fill_value=(yI[0],yI[-1]),assume_sorted=True)

    def get_radiation(self,ls):
        ''' Get radiation intensity from Allen
            Input: numpy array of wavelengths [nm]
            Output: continuum intensity at the requested frequencies
                    [erg s^-1 cm^-2 sr^-1 Hz^-1]
            Internal: None
        '''

        # Convert to lambda [micron]
        lamb = ls*1e-3

        # Intensity
        return self.fI(lamb)

######################################################################
######################################################################
######################################################################
######################################################################

class invi_class():
    ''' Class to help with the creation of a data file for the
        inversion mode of HanleRT-TIC
    '''

######################################################################
######################################################################

    def __init__(self,ofile):
        ''' Class initialization
        '''

        # Set output file path
        self.__ofile = ofile

        # Initialize variables to empty
        self.__istk = None
        self.__ilos = None
        self.__isig = 0
        self.__idif = 0

        self.__nx = None
        self.__ny = None
        self.__nl = None
        self.__nstk = None

        self.__wave = None
        self.__nl_wave = None

        self.__los = None
        self.__nx_los = None
        self.__ny_los = None

        self.__data = None
        self.__tdata = None
        self.__nx_data = None
        self.__ny_data = None
        self.__nl_data = None
        self.__nstk_data = None

        self.__sig = None
        self.__nx_sig = None
        self.__ny_sig = None
        self.__nl_sig = None
        self.__nstk_sig = None

        self.__dif = None
        self.__nx_dif = None
        self.__ny_dif = None
        self.__nl_dif = None
        self.__nstk_dif = None

######################################################################
######################################################################

    def __sanity_check(self):
        ''' Check consistent dimensions
        '''

        # Wavelength and data
        if self.__nl_wave is not None and \
           self.__tdata is not None:
            # Check dimensions
            if self.__nl_wave != self.__nl_data:
                msg = 'The lambda dimensions for ' + \
                      'lambda axis ({0}) ' + \
                      'and data ({1}) differ'
                print(msg.format(self.__nl_wave, \
                                 self.__nl_data))
                return False

        # LOS and data
        if self.__ilos is not None and \
           self.__tdata is not None:

            # If pixel-wise mu
            if self.__ilos == 1:
                # Check dimensions
                if self.__nx_los != self.__nx_data:
                    msg = 'The x dimensions for LOS ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__nx_los, \
                                     self.__nx_data))
                    return False
                if self.__ny_los != self.__ny_data:
                    msg = 'The y dimensions for LOS ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__ny_los, \
                                     self.__ny_data))
                    return False

        # Sigma and data
        if self.__isig > 0 and self.__tdata is not None:

            # Check dimensions
            if self.__nstk_sig != self.__nstk_data:
                msg = 'The Stokes dimensions for sigma ({0}) ' + \
                      'and data ({1}) differ'
                print(msg.format(self.__nstk_sig, \
                                 self.__nstk_data))
                return False
            # Pixel-wise sigma
            if self.__isig == 3 or self.__isig == 4:
                # Check dimensions
                if self.__nx_sig != self.__nx_data:
                    msg = 'The x dimensions for sigma ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__nx_sig, \
                                     self.__nx_data))
                    return False
                if self.__ny_sig != self.__ny_data:
                    msg = 'The y dimensions for sigma ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__ny_sig, \
                                     self.__ny_data))
                    return False
            # Wavelength dependent
            if self.__isig == 2 or self.__isig == 4:
                # Check dimensions
                if self.__nl_sig != self.__nl_data:
                    msg = 'The lambda dimensions for sigma ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__nl_sig, \
                                     self.__nl_data))
                    return False

        # Diffuse light and data
        if self.__idif > 0 and self.__tdata is not None:

            # Check dimensions
            if self.__nl_dif != self.__nl_data:
                msg = 'The lambda dimensions for ' + \
                      'diffuse light ({0}) ' + \
                      'and data ({1}) differ'
                print(msg.format(self.__nl_dif, \
                                 self.__nl_data))
                return False
            # Pixel-wise diffuse light
            if self.__idif == 3 or self.__idif == 4:
                # Check dimensions
                if self.__nx_dif != self.__nx_data:
                    msg = 'The x dimensions for ' + \
                          'diffuse light ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__nx_dif, \
                                     self.__nx_data))
                    return False
                if self.__ny_dif != self.__ny_data:
                    msg = 'The y dimensions for ' + \
                          'diffuse light ({0}) ' + \
                          'and data ({1}) differ'
                    print(msg.format(self.__ny_dif, \
                                     self.__ny_data))
                    return False

        return True

######################################################################
######################################################################

    def load_wavelength(self,data):
        ''' Load the wavelength in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The wavelength must be a numpy array')
            return False

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims != 1:
            print('Unexpected number of dimensions in wavelength. '+ \
                  'Got {0} instead of 1'.format(dims))
            return False

        # Define dimensions
        self.__nl_wave = data.size

        # Sanity check
        if not self.__sanity_check():
            self.__nl_wave = None
            return False

        # Save in wavelength
        self.__wave = data.copy()

        # Success
        return True

######################################################################
######################################################################

    def load_los(self,data,ix=None,iy=None):
        ''' Load the LOS in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The LOS data must be a numpy array')
            return False
        # Sanity check
        if ix is not None:
            if not isinstance(ix, int):
                print('ix must be an integer')
                return False
            jx = ix
        if iy is not None:
            if not isinstance(iy, int):
                print('iy must be an integer')
                return False
            jy = iy

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims != 1 and dims != 3:
            print('Unexpected number of dimensions in LOS. ' + \
                  'Got {0} instead of 1 or 3'.format(dims))
            return False

        #
        # Now treat the data differently depending on the
        # dimensions

        # Only one dimensions
        if dims == 1:

            # Verbose
            print('Provided a single LOS')

            # Check dims
            if data.size != 2:
                msg = 'The provided single LOS array has size ' + \
                      '{0}, but must be 2'
                print(msg.format(data.size))
                return False

            # Define dimensions
            self.__ilos = 0
            self.__nx_los = 1
            self.__ny_los = 1

            # Sanity check
            if not self.__sanity_check():
                self.__nx_los = None
                self.__ny_los = None
                self.__ilos = None
                return False

            # Save in data
            self.__los = np.empty((self.__nx_los, \
                                   self.__ny_los, 2))
            self.__los[0,0,:] = data

        # Three dimensions
        elif dims == 3:

            # Verbose
            print('Provided pixel-wise LOS')

            # Check there is info
            if ix is None and iy is None:
                msg = 'Need to provide indexes for both x and ' + \
                      'y dimensions'
                print(msg)
                return

            # Sanity
            if ix is not None:
                if jx < 0 or jx >= 3:
                    print('ix index out of bounds [0,2]')
                    return False
            if iy is not None:
                if jy < 0 or jy >= 3:
                    print('iy index out of bounds [0,2]')
                    return False
            if ix is not None and iy is not None:
                if jx == jy:
                    print('ix and iy cannot be equal')
                    return False

            # Undefined ilos
            jlos = None
            for k in range(3):
                if k == jx: continue
                if k == jy: continue
                jlos = k
                break

            # Check indexes
            if jlos is None:
                msg = 'Had a problem finding the non-spatial ' + \
                      'dimension in the LOS array'
                print(msg)
                return False
            if jx < 0 or jx > 2:
                msg = 'Got wrong x index ({0}), not in [0,2]'
                print(msg.format(jstk))
                return False
            if jy < 0 or jy > 2:
                msg = 'Got wrong y index ({0}), not in [0,2]'
                print(msg.format(jstk))
                return False
            if jlos < 0 or jlos > 3:
                msg = 'Got wrong los index ({0}), not in [0,2]'
                print(msg.format(jl))
                return False
            if jx == jy:
                print('x and y indexes cannot be equal')
                return False
            if jx == jlos:
                print('x and LOS indexes cannot be equal')
                return False
            if jy == jlos:
                print('y and LOS indexes cannot be equal')
                return False
            if shape[jlos] != 2:
                msg = 'The LOS dimension has size {0}, not 2'
                print(msg.format(shape[jstk]))
                return False

            # Define dimensions
            self.__ilos = 1
            self.__nx_los = shape[jx]
            self.__ny_los = shape[jy]

            # Sanity check
            if not self.__sanity_check():
                self.__nx_los = None
                self.__ny_los = None
                self.__ilos = None
                return False

            # Save in data
            self.__los = np.transpose(data, (jx,jy,jlos))

        # Error
        else:
            print('Somehow I reached this error for data ' + \
                  'with invalid dimensions')
            return False

        # Success
        return True

######################################################################
######################################################################

    def load_data(self,data,ix=None,iy=None,il=None,istk=None):
        ''' Load the data in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The data must be a numpy array')
            return False
        # Sanity check
        if ix is not None:
            if not isinstance(ix, int):
                print('ix must be an integer')
                return False
            jx = ix
        if iy is not None:
            if not isinstance(iy, int):
                print('iy must be an integer')
                return False
            jy = iy
        if il is not None:
            if not isinstance(il, int):
                print('il must be an integer')
                return False
            jl = il
        if istk is not None:
            if not isinstance(istk, int):
                print('istk must be an integer')
                return False
            jstk = istk

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims < 1 and dims > 4:
            print('Unexpected number of dimensions in data. ' + \
                  'Got {0} instead of 1, 2, 3, or 4'.format(dims))
            return False

        #
        # Now treat the data differently depending on the
        # dimensions

        # Only one dimensions, we need to assume that it is
        # wavelength
        if dims == 1:

            # Verbose
            print('Provided a single intensity spectrum')

            # Define dimensions
            self.__nx_data = 1
            self.__ny_data = 1
            self.__nl_data = data.size
            self.__nstk_data = 1
            self.__tdata = 1

            # Sanity check
            if not self.__sanity_check():
                self.__nx_data = None
                self.__ny_data = None
                self.__nl_data = None
                self.__nstk_data = None
                self.__tdata = None
                return False

            # Save in data
            self.__data = np.empty((self.__nx_data, \
                                    self.__ny_data, \
                                    self.__nstk_data, \
                                    self.__nl_data))
            self.__data[0,0,0,:] = data

        # Only two dimensions, we need to assume that it is
        # wavelength and Stokes
        elif dims == 2:

            # Verbose
            print('Provided a single polarized spectrum')

            # If both provided
            if istk is not None and il is not None:
                # Check indexes
                if jstk < 0 or jstk > 1:
                    msg = 'Got wrong Stokes index ({0}), not 0 or 1'
                    print(msg.format(jstk))
                    return False
                if jl < 0 or jl > 1:
                    msg = 'Got wrong lambda index ({0}), not 0 or 1'
                    print(msg.format(jl))
                    return False
                if jstk == jl:
                    print('Stokes and lambda indexes cannot be equal')
                    return False

            # Provided istk?
            if istk is None and il is None:
                if shape[0] == 4:
                    jstk = 0
                    jl = 1
                else:
                    jstk = 1
                    jl = 0
            elif istk is not None:
                if jstk == 0:
                    jl = 1
                else:
                    jl = 0
            elif il is not None:
                if jl == 0:
                    jstk == 1
                else:
                    jstk == 0

            # Check sizes
            if shape[jstk] != 4:
                msg = 'The Stokes dimension ({0}) has ' + \
                      'size {1}, expected 4'
                print(msg.format(jstk,shape[jstk]))
                return False

            # Define dimensions
            self.__nx_data = 1
            self.__ny_data = 1
            self.__nl_data = shape[jl]
            self.__nstk_data = shape[jstk]
            self.__tdata = 1

            # Sanity check
            if not self.__sanity_check():
                self.__nx_data = None
                self.__ny_data = None
                self.__nl_data = None
                self.__nstk_data = None
                self.__tdata = None
                return False

            # Save in data
            self.__data = np.empty((self.__nx_data, \
                                    self.__ny_data, \
                                    self.__nstk_data, \
                                    self.__nl_data))
            if jstk == 0:
                self.__data[0,0,:,:] = data
            else:
                self.__data[0,0,:,:] = np.transpose(data)

        # Only three dimensions, we need to assume that it is
        # only intensity
        elif dims == 3:

            # Verbose
            print('Provided a full intensity data cube')

            # Check Stokes
            if istk is not None:
                msg = 'You provided a Stokes index for ' + \
                      'a cube with dimension 3 (only intensity).' + \
                      ' Stopping just in case.'
                print(msg)
                return False

            # Check there is info
            if (ix is None and iy is None) or \
               (ix is None and il is None) or \
               (iy is None and il is None):
                msg = 'At least two indexes need to be indicated'
                print(msg)
                return False

            # Sanity
            if ix is not None:
                if jx < 0 or jx >= 3:
                    print('ix index out of bounds [0,2]')
                    return False
            if iy is not None:
                if jy < 0 or jy >= 3:
                    print('iy index out of bounds [0,2]')
                    return False
            if il is not None:
                if jl < 0 or jl >= 3:
                    print('il index out of bounds [0,2]')
                    return False
            if ix is not None and iy is not None:
                if jx == jy:
                    print('ix and iy cannot be equal')
                    return False
            if ix is not None and il is not None:
                if jx == jl:
                    print('ix and il cannot be equal')
                    return False
            if iy is not None and il is not None:
                if jy == jl:
                    print('iy and il cannot be equal')
                    return False

            # Undefined ix
            if ix is None:
                for k in range(3):
                    if k == jy: continue
                    if k == jl: continue
                    jx = k
                    break
            # Undefined iy
            if iy is None:
                for k in range(3):
                    if k == jx: continue
                    if k == jl: continue
                    jy = k
                    break
            # Undefined il
            if il is None:
                for k in range(3):
                    if k == jx: continue
                    if k == jy: continue
                    jl = k
                    break

            # Check indexes
            if jx < 0 or jx > 2:
                msg = 'Got wrong x index ({0}), not in [0,2]'
                print(msg.format(jstk))
                return False
            if jy < 0 or jy > 2:
                msg = 'Got wrong y index ({0}), not in [0,2]'
                print(msg.format(jstk))
                return False
            if jl < 0 or jl > 2:
                msg = 'Got wrong lambda index ({0}), not in [0,2]'
                print(msg.format(jl))
                return False
            if jx == jy:
                print('x and y indexes cannot be equal')
                return False
            if jx == jl:
                print('x and lambda indexes cannot be equal')
                return False
            if jy == jl:
                print('y and lambda indexes cannot be equal')
                return False

            # Define dimensions
            self.__nx_data = shape[jx]
            self.__ny_data = shape[jy]
            self.__nl_data = shape[jl]
            self.__nstk_data = 1
            self.__tdata = 1

            # Sanity check
            if not self.__sanity_check():
                self.__nx_data = None
                self.__ny_data = None
                self.__nl_data = None
                self.__nstk_data = None
                self.__tdata = None
                return False

            # Save in data
            self.__data = np.empty((self.__nx_data, \
                                    self.__ny_data, \
                                    self.__nstk_data, \
                                    self.__nl_data))
            self.__data[:,:,0,:] = np.transpose(data, (jx,jy,jl))

        # Four dimensions
        elif dims == 4:

            # Verbose
            print('Provided a full data cube')

            # Check there is info
            if (ix is None and iy is None and il is None) or \
               (ix is None and iy is None and istk is None) or \
               (ix is None and il is None and istk is None) or \
               (iy is None and il is None and istk is None):
                msg = 'At least three indexes need to be indicated'
                print(msg)
                return

            # Sanity
            if ix is not None:
                if jx < 0 or jx >= 4:
                    print('ix index out of bounds [0,3]')
                    return False
            if iy is not None:
                if jy < 0 or jy >= 4:
                    print('iy index out of bounds [0,3]')
                    return False
            if il is not None:
                if jl < 0 or jl >= 4:
                    print('il index out of bounds [0,3]')
                    return False
            if istk is not None:
                if jstk < 0 or jstk >= 4:
                    print('istk index out of bounds [0,3]')
                    return False
            if ix is not None and iy is not None:
                if jx == jy:
                    print('ix and iy cannot be equal')
                    return False
            if ix is not None and il is not None:
                if jx == jl:
                    print('ix and il cannot be equal')
                    return False
            if ix is not None and istk is not None:
                if jx == jstk:
                    print('ix and istk cannot be equal')
                    return False
            if iy is not None and il is not None:
                if jy == jl:
                    print('iy and il cannot be equal')
                    return False
            if iy is not None and istk is not None:
                if jy == jstk:
                    print('iy and istk cannot be equal')
                    return False
            if il is not None and istk is not None:
                if jl == jstk:
                    print('il and istk cannot be equal')
                    return False

            # Undefined ix
            if ix is None:
                for k in range(4):
                    if k == jy: continue
                    if k == jl: continue
                    if k == jstk: continue
                    jx = k
                    break
            # Undefined iy
            if iy is None:
                for k in range(4):
                    if k == jx: continue
                    if k == jl: continue
                    if k == jstk: continue
                    jy = k
                    break
            # Undefined il
            if il is None:
                for k in range(4):
                    if k == jx: continue
                    if k == jy: continue
                    if k == jstk: continue
                    jl = k
                    break
            # Undefined istk
            if istk is None:
                for k in range(4):
                    if k == jx: continue
                    if k == jy: continue
                    if k == jl: continue
                    jstk = k
                    break

            # Check indexes
            if jx < 0 or jx > 3:
                msg = 'Got wrong x index ({0}), not in [0,3]'
                print(msg.format(jstk))
                return False
            if jy < 0 or jy > 3:
                msg = 'Got wrong y index ({0}), not in [0,3]'
                print(msg.format(jstk))
                return False
            if jl < 0 or jl > 3:
                msg = 'Got wrong lambda index ({0}), not in [0,3]'
                print(msg.format(jl))
                return False
            if jstk < 0 or jstk > 3:
                msg = 'Got wrong Stokes index ({0}), not in [0,3]'
                print(msg.format(jl))
                return False
            if jx == jy:
                print('x and y indexes cannot be equal')
                return False
            if jx == jl:
                print('x and lambda indexes cannot be equal')
                return False
            if jx == jstk:
                print('x and Stokes indexes cannot be equal')
                return False
            if jy == jl:
                print('y and lambda indexes cannot be equal')
                return False
            if jy == jstk:
                print('y and Stokes indexes cannot be equal')
                return False
            if jl == jstk:
                print('lambda and Stokes indexes cannot be equal')
                return False
            if shape[jstk] != 4:
                msg = 'The Stokes dimension has size {0}, not 4'
                print(msg.format(shape[jstk]))
                return False

            # Define dimensions
            self.__nx_data = shape[jx]
            self.__ny_data = shape[jy]
            self.__nl_data = shape[jl]
            self.__nstk_data = shape[jstk]
            self.__tdata = 1

            # Sanity check
            if not self.__sanity_check():
                self.__nx_data = None
                self.__ny_data = None
                self.__nl_data = None
                self.__nstk_data = None
                self.__tdata = None
                return False

            # Save in data
            self.__data = np.transpose(data, (jx,jy,jstk,jl))

        # Error
        else:
            print('Somehow I reached this error for data ' + \
                  'with invalid dimensions')
            return False

        # Success
        return True

######################################################################
######################################################################

    def load_sig(self,data,ix=None,iy=None,il=None,istk=None, \
                 polarization=None):
        ''' Load the sigma in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The data must be a numpy array')
            return False
        # Sanity check
        if ix is not None:
            if not isinstance(ix, int):
                print('ix must be an integer')
                return False
            jx = ix
        if iy is not None:
            if not isinstance(iy, int):
                print('iy must be an integer')
                return False
            jy = iy
        if il is not None:
            if not isinstance(il, int):
                print('il must be an integer')
                return False
            jl = il
        if istk is not None:
            if not isinstance(istk, int):
                print('istk must be an integer')
                return False
            jstk = istk
        if polarization is not None:
            if not isinstance(polarization, bool):
                print('polarization must be an bool')
                return False
            ipol = polarization
        else:
            print('polarization is a required keyword in load_sig')
            return False

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims < 1 and dims > 4:
            print('Unexpected number of dimensions in sigma. ' + \
                  'Got {0} instead of 1, 2, 3, or 4'.format(dims))
            return False

        #
        # Polarized data
        if ipol:

            #
            # Now treat the data differently depending on the
            # dimensions

            # Only one dimensions, we need to assume that it is
            # Stokes
            if dims == 1:

                # Verbose
                print('Provided a single Stokes sigma')

                # Sanity
                if data.shape != 4:
                    msg = 'Got size {0} for Stokes, must be 4'
                    print(msg.format(data.shape))
                    return False

                # Define dimensions
                self.__isig = 1
                self.__nx_sig = 1
                self.__ny_sig = 1
                self.__nl_sig = 1
                self.__nstk_sig = 4

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in data
                self.__sig = data.copy()

            # Only two dimensions, we need to assume that it is
            # wavelength and Stokes
            elif dims == 2:

                # Verbose
                print('Provided a single polarized sigma')

                # If both provided
                if istk is not None and il is not None:
                    # Check indexes
                    if jstk < 0 or jstk > 1:
                        msg = 'Got wrong Stokes index ({0}), ' + \
                              'not 0 or 1'
                        print(msg.format(jstk))
                        return False
                    if jl < 0 or jl > 1:
                        msg = 'Got wrong lambda index ({0}), ' + \
                              'not 0 or 1'
                        print(msg.format(jl))
                        return False
                    if jstk == jl:
                        print('Stokes and lambda indexes ' + \
                              'cannot be equal')
                        return False

                # Provided istk?
                if istk is None and il is None:
                    if shape[0] == 4:
                        jstk = 0
                        jl = 1
                    else:
                        jstk = 1
                        jl = 0
                elif istk is not None:
                    if jstk == 0:
                        jl = 1
                    else:
                        jl = 0
                elif il is not None:
                    if jl == 0:
                        jstk == 1
                    else:
                        jstk == 0

                # Check sizes
                if shape[jstk] != 4:
                    msg = 'The Stokes dimension ({0}) has ' + \
                          'size {1}, expected 4'
                    print(msg.format(jstk,shape[jstk]))
                    return False

                # Define dimensions
                self.__isig = 2
                self.__nx_sig = 1
                self.__ny_sig = 1
                self.__nl_sig = shape[jl]
                self.__nstk_sig = shape[jstk]

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in sigma
                self.__sig = np.transpose(data, (jstk,jl))

            # Only three dimensions, we need to assume that it is
            # wavelength independent
            elif dims == 3:

                # Verbose
                print('Provided a pixel dependent constant sigma')

                # Check there is info
                if (ix is None and iy is None) or \
                   (ix is None and istk is None) or \
                   (iy is None and istk is None):
                    msg = 'At least two indexes need to be indicated'
                    print(msg)
                    return False

                # Sanity
                if ix is not None:
                    if jx < 0 or jx >= 3:
                        print('ix index out of bounds [0,2]')
                        return False
                if iy is not None:
                    if jy < 0 or jy >= 3:
                        print('iy index out of bounds [0,2]')
                        return False
                if istk is not None:
                    if jstk < 0 or jstk >= 3:
                        print('istk index out of bounds [0,2]')
                        return False
                if ix is not None and iy is not None:
                    if jx == jy:
                        print('ix and iy cannot be equal')
                        return False
                if ix is not None and istk is not None:
                    if jx == jstk:
                        print('ix and istk cannot be equal')
                        return False
                if iy is not None and istk is not None:
                    if jy == jstk:
                        print('iy and istk cannot be equal')
                        return False

                # Undefined ix
                if ix is None:
                    for k in range(3):
                        if k == jy: continue
                        if k == jstk: continue
                        jx = k
                        break
                # Undefined iy
                if iy is None:
                    for k in range(3):
                        if k == jx: continue
                        if k == jstk: continue
                        jy = k
                        break
                # Undefined istk
                if istk is None:
                    for k in range(3):
                        if k == jx: continue
                        if k == jy: continue
                        jstk = k
                        break

                # Check indexes
                if jx < 0 or jx > 2:
                    msg = 'Got wrong x index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jy < 0 or jy > 2:
                    msg = 'Got wrong y index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jstk < 0 or jstk > 2:
                    msg = 'Got wrong Stokes index ({0}), not in [0,2]'
                    print(msg.format(jl))
                    return False
                if jx == jy:
                    print('x and y indexes cannot be equal')
                    return False
                if jx == jstk:
                    print('x and Stokes indexes cannot be equal')
                    return False
                if jy == jstk:
                    print('y and Stokes indexes cannot be equal')
                    return False
                if shape[jstk] != 4:
                    msg = 'The Stokes dimension ({0}) has ' + \
                          'size {1}, expected 4'
                    print(msg.format(jstk,shape[jstk]))
                    return False

                # Define dimensions
                self.__isig = 3
                self.__nx_sig = shape[jx]
                self.__ny_sig = shape[jy]
                self.__nl_sig = 1
                self.__nstk_sig = shape[jstk]

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in sig
                self.__sig = np.empty((self.__nx_sig, \
                                       self.__ny_sig, \
                                       self.__nstk_sig, \
                                       self.__nl_sig))
                self.__sig[:,:,:,0] = np.transpose(data, (jx,jy,jstk))

            # Four dimensions
            elif dims == 4:

                # Verbose
                print('Provided a full sigma cube')

                # Check there is info
                if (ix is None and iy is None and il is None) or \
                   (ix is None and iy is None and istk is None) or \
                   (ix is None and il is None and istk is None) or \
                   (iy is None and il is None and istk is None):
                    msg = 'At least three indexes need to be ' + \
                          'indicated'
                    print(msg)
                    return

                # Sanity
                if ix is not None:
                    if jx < 0 or jx >= 4:
                        print('ix index out of bounds [0,3]')
                        return False
                if iy is not None:
                    if jy < 0 or jy >= 4:
                        print('iy index out of bounds [0,3]')
                        return False
                if il is not None:
                    if jl < 0 or jl >= 4:
                        print('il index out of bounds [0,3]')
                        return False
                if istk is not None:
                    if jstk < 0 or jstk >= 4:
                        print('istk index out of bounds [0,3]')
                        return False
                if ix is not None and iy is not None:
                    if jx == jy:
                        print('ix and iy cannot be equal')
                        return False
                if ix is not None and il is not None:
                    if jx == jl:
                        print('ix and il cannot be equal')
                        return False
                if ix is not None and istk is not None:
                    if jx == jstk:
                        print('ix and istk cannot be equal')
                        return False
                if iy is not None and il is not None:
                    if jy == jl:
                        print('iy and il cannot be equal')
                        return False
                if iy is not None and istk is not None:
                    if jy == jstk:
                        print('iy and istk cannot be equal')
                        return False
                if il is not None and istk is not None:
                    if jl == jstk:
                        print('il and istk cannot be equal')
                        return False

                # Undefined ix
                if ix is None:
                    for k in range(4):
                        if k == jy: continue
                        if k == jl: continue
                        if k == jstk: continue
                        jx = k
                        break
                # Undefined iy
                if iy is None:
                    for k in range(4):
                        if k == jx: continue
                        if k == jl: continue
                        if k == jstk: continue
                        jy = k
                        break
                # Undefined il
                if il is None:
                    for k in range(4):
                        if k == jx: continue
                        if k == jy: continue
                        if k == jstk: continue
                        jl = k
                        break
                # Undefined istk
                if istk is None:
                    for k in range(4):
                        if k == jx: continue
                        if k == jy: continue
                        if k == jl: continue
                        jstk = k
                        break

                # Check indexes
                if jx < 0 or jx > 3:
                    msg = 'Got wrong x index ({0}), not in [0,3]'
                    print(msg.format(jstk))
                    return False
                if jy < 0 or jy > 3:
                    msg = 'Got wrong y index ({0}), not in [0,3]'
                    print(msg.format(jstk))
                    return False
                if jl < 0 or jl > 3:
                    msg = 'Got wrong lambda index ({0}), not in [0,3]'
                    print(msg.format(jl))
                    return False
                if jstk < 0 or jstk > 3:
                    msg = 'Got wrong Stokes index ({0}), not in [0,3]'
                    print(msg.format(jl))
                    return False
                if jx == jy:
                    print('x and y indexes cannot be equal')
                    return False
                if jx == jl:
                    print('x and lambda indexes cannot be equal')
                    return False
                if jx == jstk:
                    print('x and Stokes indexes cannot be equal')
                    return False
                if jy == jl:
                    print('y and lambda indexes cannot be equal')
                    return False
                if jy == jstk:
                    print('y and Stokes indexes cannot be equal')
                    return False
                if jl == jstk:
                    print('lambda and Stokes indexes cannot be equal')
                    return False
                if shape[jstk] != 4:
                    msg = 'The Stokes dimension has size {0}, not 4'
                    print(msg.format(shape[jstk]))
                    return False

                # Define dimensions
                self.__isig = 4
                self.__nx_sig = shape[jx]
                self.__ny_sig = shape[jy]
                self.__nl_sig = shape[jl]
                self.__nstk_sig = shape[jstk]

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in sigma
                self.__sig = np.transpose(data, (jx,jy,jstk,jl))

            # Error
            else:
                print('Somehow I reached this error for sigma ' + \
                      'with invalid dimensions')
                return False
        #
        # Unpolarized
        else:

            #
            # Now treat the data differently depending on the
            # dimensions

            # Only one dimensions, we need to assume that it is
            # a number
            if dims == 1:

                # Single number
                if data.size == 1:

                    # Verbose
                    print('Provided a single intensity sigma')

                    # Define dimensions
                    self.__isig = 1
                    self.__nx_sig = 1
                    self.__ny_sig = 1
                    self.__nl_sig = 1
                    self.__nstk_sig = 1

                else:

                    # Verbose
                    print('Provided wavelength dependent ' + \
                          'intensity sigma')

                    # Define dimensions
                    self.__isig = 2
                    self.__nx_sig = 1
                    self.__ny_sig = 1
                    self.__nl_sig = data.size
                    self.__nstk_sig = 1

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in data
                self.__sig = data.copy()

            # Only two dimensions, we need to assume that it is
            # wavelength independent
            elif dims == 2:

                # Verbose
                print('Provided a pixel dependent ' + \
                      'constant intensity sigma')

                # Check there is info
                if ix is None and iy is None:
                    msg = 'At least one index needs to be indicated'
                    print(msg)
                    return False

                # Sanity
                if ix is not None:
                    if jx < 0 or jx >= 3:
                        print('ix index out of bounds [0,2]')
                        return False
                if iy is not None:
                    if jy < 0 or jy >= 3:
                        print('iy index out of bounds [0,2]')
                        return False
                if ix is not None and iy is not None:
                    if jx == jy:
                        print('ix and iy cannot be equal')
                        return False

                # Undefined ix
                if ix is None:
                    for k in range(3):
                        if k == jy: continue
                        jx = k
                        break
                # Undefined iy
                if iy is None:
                    for k in range(3):
                        if k == jx: continue
                        jy = k
                        break

                # Check indexes
                if jx < 0 or jx > 2:
                    msg = 'Got wrong x index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jy < 0 or jy > 2:
                    msg = 'Got wrong y index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jx == jy:
                    print('x and y indexes cannot be equal')
                    return False

                # Define dimensions
                self.__isig = 3
                self.__nx_sig = shape[jx]
                self.__ny_sig = shape[jy]
                self.__nl_sig = 1
                self.__nstk_sig = 1

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in sig
                self.__sig = np.empty((self.__nx_sig, \
                                       self.__ny_sig, \
                                       self.__nstk_sig, \
                                       self.__nl_sig))
                self.__sig[:,:,0,0] = np.transpose(data, (jx,jy))

            # Three dimensions
            elif dims == 3:

                # Verbose
                print('Provided a full intensity sigma cube')

                # Check there is info
                if (ix is None and iy is None) or \
                   (ix is None and il is None) or \
                   (iy is None and il is None):
                    msg = 'At least two indexes need to be ' + \
                          'indicated'
                    print(msg)
                    return

                # Sanity
                if ix is not None:
                    if jx < 0 or jx >= 3:
                        print('ix index out of bounds [0,2]')
                        return False
                if iy is not None:
                    if jy < 0 or jy >= 3:
                        print('iy index out of bounds [0,2]')
                        return False
                if il is not None:
                    if jl < 0 or jl >= 3:
                        print('il index out of bounds [0,2]')
                        return False
                if ix is not None and iy is not None:
                    if jx == jy:
                        print('ix and iy cannot be equal')
                        return False
                if ix is not None and il is not None:
                    if jx == jl:
                        print('ix and il cannot be equal')
                        return False
                if iy is not None and il is not None:
                    if jy == jl:
                        print('iy and il cannot be equal')
                        return False

                # Undefined ix
                if ix is None:
                    for k in range(3):
                        if k == jy: continue
                        if k == jl: continue
                        jx = k
                        break
                # Undefined iy
                if iy is None:
                    for k in range(3):
                        if k == jx: continue
                        if k == jl: continue
                        jy = k
                        break
                # Undefined il
                if il is None:
                    for k in range(3):
                        if k == jx: continue
                        if k == jy: continue
                        jl = k
                        break

                # Check indexes
                if jx < 0 or jx >= 3:
                    msg = 'Got wrong x index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jy < 0 or jy >= 3:
                    msg = 'Got wrong y index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jl < 0 or jl >= 3:
                    msg = 'Got wrong lambda index ({0}), not in [0,2]'
                    print(msg.format(jl))
                    return False
                if jx == jy:
                    print('x and y indexes cannot be equal')
                    return False
                if jx == jl:
                    print('x and lambda indexes cannot be equal')
                    return False
                if jy == jl:
                    print('y and lambda indexes cannot be equal')
                    return False

                # Define dimensions
                self.__isig = 4
                self.__nx_sig = shape[jx]
                self.__ny_sig = shape[jy]
                self.__nl_sig = shape[jl]
                self.__nstk_sig = 1

                # Sanity check
                if not self.__sanity_check():
                    self.__isig = 0
                    self.__nx_sig = None
                    self.__ny_sig = None
                    self.__nl_sig = None
                    self.__nstk_sig = None
                    return False

                # Save in sigma
                self.__sig = np.empty((self.__nx_sig, \
                                       self.__ny_sig, \
                                       self.__nstk_sig, \
                                       self.__nl_sig))
                self.__sig[:,:,0,:] = np.transpose(data, \
                                                   (jx,jy,jl))

            # Error
            else:
                print('Somehow I reached this error for sigma ' + \
                      'with invalid dimensions')
                return False

        # Success
        return True

######################################################################
######################################################################

    def load_dif(self,data,ix=None,iy=None,il=None,istk=None, \
                 polarization=None):
        ''' Load the diffuse light in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The data must be a numpy array')
            return False
        # Sanity check
        if ix is not None:
            if not isinstance(ix, int):
                print('ix must be an integer')
                return False
            jx = ix
        if iy is not None:
            if not isinstance(iy, int):
                print('iy must be an integer')
                return False
            jy = iy
        if il is not None:
            if not isinstance(il, int):
                print('il must be an integer')
                return False
            jl = il
        if istk is not None:
            if not isinstance(istk, int):
                print('istk must be an integer')
                return False
            jstk = istk
        if polarization is not None:
            if not isinstance(polarization, bool):
                print('polarization must be an bool')
                return False
            ipol = polarization
        else:
            print('polarization is a required keyword in load_sig')
            return False

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims < 1 and dims > 4:
            print('Unexpected number of dimensions in sigma. ' + \
                  'Got {0} instead of 1, 2, 3, or 4'.format(dims))
            return False

        #
        # Polarized data
        if ipol:

            #
            # Now treat the data differently depending on the
            # dimensions

            # Only two dimensions, we need to assume that it is
            # wavelength and Stokes
            if dims == 2:

                # Verbose
                print('Provided a single polarized ' + \
                      'diffuse light profile')

                # If both provided
                if istk is not None and il is not None:
                    # Check indexes
                    if jstk < 0 or jstk > 1:
                        msg = 'Got wrong Stokes index ({0}), ' + \
                              'not 0 or 1'
                        print(msg.format(jstk))
                        return False
                    if jl < 0 or jl > 1:
                        msg = 'Got wrong lambda index ({0}), ' + \
                              'not 0 or 1'
                        print(msg.format(jl))
                        return False
                    if jstk == jl:
                        print('Stokes and lambda indexes ' + \
                              'cannot be equal')
                        return False

                # Provided istk?
                if istk is None and il is None:
                    if shape[0] == 4:
                        jstk = 0
                        jl = 1
                    else:
                        jstk = 1
                        jl = 0
                elif istk is not None:
                    if jstk == 0:
                        jl = 1
                    else:
                        jl = 0
                elif il is not None:
                    if jl == 0:
                        jstk == 1
                    else:
                        jstk == 0

                # Check sizes
                if shape[jstk] != 4:
                    msg = 'The Stokes dimension ({0}) has ' + \
                          'size {1}, expected 4'
                    print(msg.format(jstk,shape[jstk]))
                    return False

                # Define dimensions
                self.__idif = 2
                self.__nx_dif = 1
                self.__ny_dif = 1
                self.__nl_dif = shape[jl]
                self.__nstk_dif = shape[jstk]

                # Sanity check
                if not self.__sanity_check():
                    self.__idif = 0
                    self.__nx_dif = None
                    self.__ny_dif = None
                    self.__nl_dif = None
                    self.__nstk_dif = None
                    return False

                # Save in diffuse light
                self.__dif = np.transpose(data, (jstk,jl))

            # Four dimensions
            elif dims == 4:

                # Verbose
                print('Provided a full Stokes diffuse light cube')

                # Check there is info
                if (ix is None and iy is None and il is None) or \
                   (ix is None and iy is None and istk is None) or \
                   (ix is None and il is None and istk is None) or \
                   (iy is None and il is None and istk is None):
                    msg = 'At least three indexes need to be ' + \
                          'indicated'
                    print(msg)
                    return

                # Sanity
                if ix is not None:
                    if jx < 0 or jx >= 4:
                        print('ix index out of bounds [0,3]')
                        return False
                if iy is not None:
                    if jy < 0 or jy >= 4:
                        print('iy index out of bounds [0,3]')
                        return False
                if il is not None:
                    if jl < 0 or jl >= 4:
                        print('il index out of bounds [0,3]')
                        return False
                if istk is not None:
                    if jstk < 0 or jstk >= 4:
                        print('istk index out of bounds [0,3]')
                        return False
                if ix is not None and iy is not None:
                    if jx == jy:
                        print('ix and iy cannot be equal')
                        return False
                if ix is not None and il is not None:
                    if jx == jl:
                        print('ix and il cannot be equal')
                        return False
                if ix is not None and istk is not None:
                    if jx == jstk:
                        print('ix and istk cannot be equal')
                        return False
                if iy is not None and il is not None:
                    if jy == jl:
                        print('iy and il cannot be equal')
                        return False
                if iy is not None and istk is not None:
                    if jy == jstk:
                        print('iy and istk cannot be equal')
                        return False
                if il is not None and istk is not None:
                    if jl == jstk:
                        print('il and istk cannot be equal')
                        return False

                # Undefined ix
                if ix is None:
                    for k in range(4):
                        if k == jy: continue
                        if k == jl: continue
                        if k == jstk: continue
                        jx = k
                        break
                # Undefined iy
                if iy is None:
                    for k in range(4):
                        if k == jx: continue
                        if k == jl: continue
                        if k == jstk: continue
                        jy = k
                        break
                # Undefined il
                if il is None:
                    for k in range(4):
                        if k == jx: continue
                        if k == jy: continue
                        if k == jstk: continue
                        jl = k
                        break
                # Undefined istk
                if istk is None:
                    for k in range(4):
                        if k == jx: continue
                        if k == jy: continue
                        if k == jl: continue
                        jstk = k
                        break

                # Check indexes
                if jx < 0 or jx > 3:
                    msg = 'Got wrong x index ({0}), not in [0,3]'
                    print(msg.format(jstk))
                    return False
                if jy < 0 or jy > 3:
                    msg = 'Got wrong y index ({0}), not in [0,3]'
                    print(msg.format(jstk))
                    return False
                if jl < 0 or jl > 3:
                    msg = 'Got wrong lambda index ({0}), not in [0,3]'
                    print(msg.format(jl))
                    return False
                if jstk < 0 or jstk > 3:
                    msg = 'Got wrong Stokes index ({0}), not in [0,3]'
                    print(msg.format(jl))
                    return False
                if jx == jy:
                    print('x and y indexes cannot be equal')
                    return False
                if jx == jl:
                    print('x and lambda indexes cannot be equal')
                    return False
                if jx == jstk:
                    print('x and Stokes indexes cannot be equal')
                    return False
                if jy == jl:
                    print('y and lambda indexes cannot be equal')
                    return False
                if jy == jstk:
                    print('y and Stokes indexes cannot be equal')
                    return False
                if jl == jstk:
                    print('lambda and Stokes indexes cannot be equal')
                    return False
                if shape[jstk] != 4:
                    msg = 'The Stokes dimension has size {0}, not 4'
                    print(msg.format(shape[jstk]))
                    return False

                # Define dimensions
                self.__idif = 4
                self.__nx_dif = shape[jx]
                self.__ny_dif = shape[jy]
                self.__nl_dif = shape[jl]
                self.__nstk_dif = shape[jstk]

                # Sanity check
                if not self.__sanity_check():
                    self.__idif = 0
                    self.__nx_dif = None
                    self.__ny_dif = None
                    self.__nl_dif = None
                    self.__nstk_dif = None
                    return False

                # Save in diffuse light
                self.__dif = np.transpose(data, (jx,jy,jstk,jl))

            # Error
            else:
                print('Somehow I reached this error for diffuse ' + \
                      'light with invalid dimensions')
                return False
        #
        # Unpolarized
        else:

            #
            # Now treat the data differently depending on the
            # dimensions

            # Only one dimensions, we need to assume that it is
            # a number
            if dims == 1:

                # Verbose
                print('Provided single intensity diffuse light ' + \
                      'profile')

                # Define dimensions
                self.__idif = 1
                self.__nx_dif = 1
                self.__ny_dif = 1
                self.__nl_dif = data.size
                self.__nstk_dif = 1

                # Sanity check
                if not self.__sanity_check():
                    self.__idif = 0
                    self.__nx_dif = None
                    self.__ny_dif = None
                    self.__nl_dif = None
                    self.__nstk_dif = None
                    return False

                # Save in data
                self.__dif = data.copy()

            # Three dimensions
            elif dims == 3:

                # Verbose
                print('Provided a full intensity diffuse ' + \
                      'light profile cube')

                # Check there is info
                if (ix is None and iy is None) or \
                   (ix is None and il is None) or \
                   (iy is None and il is None):
                    msg = 'At least two indexes need to be ' + \
                          'indicated'
                    print(msg)
                    return

                # Sanity
                if ix is not None:
                    if jx < 0 or jx >= 3:
                        print('ix index out of bounds [0,2]')
                        return False
                if iy is not None:
                    if jy < 0 or jy >= 3:
                        print('iy index out of bounds [0,2]')
                        return False
                if il is not None:
                    if jl < 0 or jl >= 3:
                        print('il index out of bounds [0,2]')
                        return False
                if ix is not None and iy is not None:
                    if jx == jy:
                        print('ix and iy cannot be equal')
                        return False
                if ix is not None and il is not None:
                    if jx == jl:
                        print('ix and il cannot be equal')
                        return False
                if iy is not None and il is not None:
                    if jy == jl:
                        print('iy and il cannot be equal')
                        return False

                # Undefined ix
                if ix is None:
                    for k in range(3):
                        if k == jy: continue
                        if k == jl: continue
                        jx = k
                        break
                # Undefined iy
                if iy is None:
                    for k in range(3):
                        if k == jx: continue
                        if k == jl: continue
                        jy = k
                        break
                # Undefined il
                if il is None:
                    for k in range(3):
                        if k == jx: continue
                        if k == jy: continue
                        jl = k
                        break

                # Check indexes
                if jx < 0 or jx >= 3:
                    msg = 'Got wrong x index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jy < 0 or jy >= 3:
                    msg = 'Got wrong y index ({0}), not in [0,2]'
                    print(msg.format(jstk))
                    return False
                if jl < 0 or jl >= 3:
                    msg = 'Got wrong lambda index ({0}), not in [0,2]'
                    print(msg.format(jl))
                    return False
                if jx == jy:
                    print('x and y indexes cannot be equal')
                    return False
                if jx == jl:
                    print('x and lambda indexes cannot be equal')
                    return False
                if jy == jl:
                    print('y and lambda indexes cannot be equal')
                    return False

                # Define dimensions
                self.__idif = 3
                self.__nx_dif = shape[jx]
                self.__ny_dif = shape[jy]
                self.__nl_dif = shape[jl]
                self.__nstk_dif = 1

                # Sanity check
                if not self.__sanity_check():
                    self.__idif = 0
                    self.__nx_dif = None
                    self.__ny_dif = None
                    self.__nl_dif = None
                    self.__nstk_dif = None
                    return False

                # Save in diffuse light
                self.__dif = np.empty((self.__nx_dif, \
                                       self.__ny_dif, \
                                       self.__nstk_dif, \
                                       self.__nl_dif))
                self.__dif[:,:,0,:] = np.transpose(data, \
                                                   (jx,jy,jl))

            # Error
            else:
                print('Somehow I reached this error for diffuse ' + \
                      'light profile with invalid dimensions')
                return False

        # Success
        return True

######################################################################
######################################################################

    def create_file(self):
        ''' Create a file with the available data
        '''

        # A last sanity check
        if not self.__sanity_check():
            return False

        # Check there is LOS
        if self.__ilos is None:
            print('It is required to specify the line of sight (LOS)')
            return False

        # Check there is wavelength axis
        if self.__wave is None:
            print('There is no wavelength loaded!')
            return False

        # Check there is data
        if self.__data is None:
            print('There is no data loaded!')
            return False

        # Check if stokes
        if self.__nstk_data > 1:
            istk = 1
        else:
            istk = 0
        nstk = istk*3

        # Simplify dimensions
        self.__nx = self.__nx_data
        self.__ny = self.__ny_data
        self.__nl = self.__nl_wave
        self.__nstk = self.__nstk_data

        #
        # Try to create file
        try:

            # Open file
            f = open(self.__ofile,'wb')

            # Label
            f.write('invi'.encode())

            # Dimensions
            f.write(struct.pack('i',self.__nx))
            f.write(struct.pack('i',self.__ny))
            f.write(struct.pack('i',self.__nl))

            # Info Stokes
            f.write(struct.pack('i',istk))

            # Info mu
            f.write(struct.pack('i',self.__ilos))

            # Info sigma
            f.write(struct.pack('i',self.__isig))

            # Info diffuse light
            if self.__idif == 2 and istk == 0:
                f.write(struct.pack('i',1))
            elif self.__idif == 4 and istk == 0:
                f.write(struct.pack('i',3))
            else:
                f.write(struct.pack('i',self.__idif))

            # Wavelengths
            f.write(struct.pack('d'*self.__nl, *self.__wave))

            # Constant LOS?
            if self.__ilos == 0:
                f.write(struct.pack('dd', *(self.__los.ravel())))

            # Unpolarized
            if istk == 0:

                # Constant sigma
                if self.__isig == 1:

                    # Write
                    f.write(struct.pack('d',self.__sig))

                # Wavelength dependent sigma
                elif self.__isig == 2:

                    # Write
                    f.write(struct.pack('d'*self.__nl, \
                                        *self.__sig))

            # Polarized
            else:

                # Constant
                if self.__isig == 1:

                    # Write
                    f.write(struct.pack('dddd',*self.__sig))

                # Wavelength dependent sigma
                elif self.__isig == 2:

                    # Write
                    f.write(struct.pack('d'*self.__nl*4, \
                                        *self.__sig.flatten()))

            # Constant intensity diffuse light
            if self.__idif == 1:

                # Write
                f.write(struct.pack('d'*self.__nl, *self.__dif))

            # Constant polarization diffuse light
            elif self.__idif == 2:

                # If unpolarized
                if istk == 0:

                    # Write
                    f.write(struct.pack('d'*self.__nl, \
                                        *self.__dif[0,:]))

                # Polarized
                else:

                    # Write
                    f.write(struct.pack('d'*self.__nl*4, \
                                        *self.__dif.flatten()))

            # For each pixel
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Write LOS
                    if self.__ilos == 1:
                        f.write(struct.pack('dd', \
                                            *self.__los[ix,iy,:]))

                    # Write data
                    f.write(struct.pack('d'*self.__nl*self.__nstk, \
                                   *self.__data[ix,iy,:,:].flatten()))

                    # If sigma
                    if self.__isig > 2:

                        # Write
                        f.write(struct.pack('d'*self.__nl_sig* \
                                                self.__nstk_sig, \
                                    *self.__sig[ix,iy,:,:].flatten()))

                    # If diffuse light
                    if self.__idif > 2:

                        # Unpolarized data
                        if istk == 0:

                            # Write
                            f.write(struct.pack('d'*self.__nl_dif, \
                                              *self.__dif[ix,iy,0,:]))

                        # Polarized data
                        else:

                            # Write
                            f.write(struct.pack('d'*self.__nl_dif* \
                                                self.__nstk_dif, \
                                    *self.__dif[ix,iy,:,:].flatten()))

            # Close file
            f.close()

        except:
            raise

        # Success
        return True

######################################################################
######################################################################
######################################################################
######################################################################

class atmo15D_class():
    ''' Class to help with the creation of an atmospheric model for
        1.5D synthesis in HanleRT-TIC
    '''

######################################################################
######################################################################

    def __init__(self,ofile):
        ''' Class initialization
        '''

        # Set output file path
        self.__ofile = ofile

        # Initialize variables to empty
        self.__data = []
        for i in range(24):
            self.__data.append(None)

        self.__nx = None
        self.__ny = None
        self.__nz = None
        self.__prec = None
        self.__ltau = False

        # Index dictionary
        self.__idx_dic = {'h': 0, \
                          'tau': 1, \
                          'ltau': 1, \
                          'chi': 2, \
                          't': 3, \
                          'pg': 4, \
                          'rho': 5, \
                          'bx': 6, \
                          'by': 7, \
                          'bz': 8, \
                          'vx': 9, \
                          'vy': 10, \
                          'vz': 11, \
                          'vmi': 12, \
                          'pe': 13, \
                          'ne': 14, \
                          'nht': 15, \
                          'nha': 16, \
                          'nhm': 17, \
                          'nh0': 18, \
                          'nh1': 19, \
                          'nh2': 20, \
                          'nh3': 21, \
                          'nh4': 22, \
                          'np': 23}

        # Descriptions dictionary
        self.__des_dic = {'h': 'height [km]', \
                          'tau/ltau': 'optical depth ' + \
                                      '(not/yes logarithmic)', \
                          't': 'temperature [K]', \
                          'pg': 'Gas pressure [dyn cm^-2]', \
                          'rho': 'Density [cm^-3]',  \
                          'bx': 'Magnetic field X component [G]', \
                          'by': 'Magnetic field Y component [G]', \
                          'bz': 'Magnetic field Z component [G]', \
                          'vx': 'Velocity X component [km s^-1]', \
                          'vy': 'Velocity Y component [km s^-1]', \
                          'vz': 'Velocity Z component [km s^-1]', \
                          'vmi': 'Microturbulent velocity', \
                          'pe': 'Electron pressure [dyn cm^-2]', \
                          'ne': 'Electron number density [cm^-3]', \
                          'nht': 'Total hydrogen number ' + \
                                 'density [cm^-3]', \
                          'nha': 'Atomic hydrogen number ' + \
                                 'density [cm^-3]', \
                          'nhm': 'H- number density [cm^-3]', \
                          'nh0': 'Ground level hydrogen number ' + \
                                 'density [cm^-3]', \
                          'nh1': 'First excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'nh2': 'Second excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'nh3': 'Third excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'nh4': 'Fourth excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'np': 'Proton number density [cm^-3]'}

######################################################################
######################################################################

    def __sanity_check(self,shape=None):
        ''' Check consistent dimensions
        '''

        # Precision
        if self.__prec is not None:
            if self.__prec != np.float32 and \
               self.__prec != np.float64:
                msg = 'The precision ({0}) ' + \
                      'is not single nor double'
                print(msg.format(self.__prec))
                return False

        # If specified shape
        if shape is not None:
            if self.__nx != shape[0]:
                msg = 'The x dimensions ({0}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({1})'
                print(msg.format(shape[0],self.__nx))
                return False
            if self.__ny != shape[1]:
                msg = 'The y dimensions ({0}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({1})'
                print(msg.format(shape[1],self.__ny))
                return False
            if self.__nz != shape[2]:
                msg = 'The z dimensions ({0}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({1})'
                print(msg.format(shape[2],self.__nz))
                return False
            return True

        # For each variable
        for var,svar in zip(self.__data, \
            ['z','tau','chi','T','Pg','rho','Bx','By','Bz', \
             'vx','vy','vz','v_mi','Pe-','ne-','nH','nHa','nH-', \
             'nH_0','nH_1','nH_2','nH_3','nH_4','np+']):

            # Skip undefined
            if var is None: continue

            # Check dimensions
            if self.__nx != var.shape[0]:
                msg = 'The x dimensions for variable {0} ({1}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({2})'
                print(msg.format(svar,var.shape[0],self.__nx))
                return False
            # Check dimensions
            if self.__ny != var.shape[1]:
                msg = 'The y dimensions for variable {0} ({1}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({2})'
                print(msg.format(svar,var.shape[1],self.__ny))
                return False
            # Check dimensions
            if self.__nz != var.shape[2]:
                msg = 'The x dimensions for variable {0} ({1}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({2})'
                print(msg.format(svar,var.shape[2],self.__nz))
                return False

        # Fine
        return True

######################################################################
######################################################################

    def __check_vars(self):
        ''' Check that all necessary variables are loaded
        '''
        if self.__data[0] is None and self.__data[1] is None:
            print('Either the height or the optical depth must be ' +\
                  'defined')
            return False
        if self.__data[0] is not None:
            diff = np.min(self.__data[0][:,:,1] - \
                          self.__data[0][:,:,0])
            if diff < 0.:
                self.__orderz = True
            elif diff > 0:
                self.__orderz = False
            else:
                print('There cannot be two equal height values')
                return False
        else:
            self.__orderz = None
        if self.__data[1] is not None:
            diff = np.min(self.__data[0][:,:,1] - \
                          self.__data[0][:,:,0])
            if diff > 0.:
                self.__ordert = True
            elif diff < 0:
                self.__ordert = False
            else:
                print('There cannot be two equal optical depths')
                return False
        else:
            self.__ordert = None
        if self.__orderz is not None and self.__ordert is not None:
            if (self.__orderz and not self.__ordert) or \
               (self.__ordert and not self.__orderz):
                print('The ordering of the height and tau axes ' + \
                      'is not consistent')
                return False
        if self.__orderz is None: self.__orderz = self.__ordert
        if self.__data[3] is None:
            print('The temperature must be defined')
            return False
        if self.__data[4] is None and \
           self.__data[5] is None and \
           self.__data[13] is None and \
           self.__data[14] is None:
            print('Either the gass pressure, the density, the ' + \
                  'electron pressure, or the electron number ' + \
                  'density must be defined')
            return False
        return True

######################################################################
######################################################################

    def variable_list(self):
        ''' Return the list of available variables
        '''
        for key in self.__des_dic:
            print(f'- {key:3s}: {self.__des_dic[key]}')

######################################################################
######################################################################

    def set_precision(self,prec):
        ''' Set the kind of precision
        '''
        if not isinstance(prec,str):
            msg = 'The input to set_precision must be a string'
            print(msg)
            return False
        if prec.lower() == 's':
            self.__prec = np.float32
            return True
        elif prec.lower() == 'd':
            self.__prec = np.float64
            return True
        else:
            msg = 'The input to set_precision must be either "s" ' + \
                  'for single precision or "d" for double precision'
            return False

######################################################################
######################################################################

    def load_data(self,data,var,ix=None,iy=None,iz=None):
        ''' Load the data in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The data must be a numpy array')
            return False
        # Sanity check
        if not isinstance(var, str):
            print('Second argument must be a string')
            return False
        if ix is not None:
            if not isinstance(ix, int):
                print('ix must be an integer')
                return False
            jx = ix
        if iy is not None:
            if not isinstance(iy, int):
                print('iy must be an integer')
                return False
            jy = iy
        if iz is not None:
            if not isinstance(iz, int):
                print('iz must be an integer')
                return False
            jz = iz
        if var.lower() not in self.__idx_dic:
            print(f'Variable {var} not in the list, check with ' + \
                  'variable_list() method')
            return False

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims != 3:
            print('Unexpected number of dimensions in data. ' + \
                  'Got {0} instead of 3'.format(dims))
            return False

        # Check there is info
        if (ix is None and iy is None) or \
           (ix is None and iz is None) or \
           (iy is None and iz is None):
            msg = 'At least two indexes need to be indicated'
            print(msg)
            return False

        # Sanity
        if ix is not None:
            if jx < 0 or jx >= 3:
                print('ix index out of bounds [0,2]')
                return False
        if iy is not None:
            if jy < 0 or jy >= 3:
                print('iy index out of bounds [0,2]')
                return False
        if iz is not None:
            if jz < 0 or jz >= 3:
                print('iz index out of bounds [0,2]')
                return False
        if ix is not None and iy is not None:
            if jx == jy:
                print('ix and iy cannot be equal')
                return False
        if ix is not None and iz is not None:
            if jx == jz:
                print('ix and iz cannot be equal')
                return False
        if iy is not None and iz is not None:
            if jy == jz:
                print('iy and iz cannot be equal')
                return False

        # Undefined ix
        if ix is None:
            for k in range(3):
                if k == jy: continue
                if k == jz: continue
                jx = k
                break
        # Undefined iy
        if iy is None:
            for k in range(3):
                if k == jx: continue
                if k == jz: continue
                jy = k
                break
        # Undefined iz
        if iz is None:
            for k in range(3):
                if k == jx: continue
                if k == jy: continue
                jz = k
                break

        # Check indexes
        if jx < 0 or jx > 2:
            msg = 'Got wrong x index ({0}), not in [0,2]'
            print(msg.format(jstk))
            return False
        if jy < 0 or jy > 2:
            msg = 'Got wrong y index ({0}), not in [0,2]'
            print(msg.format(jstk))
            return False
        if jz < 0 or jz > 2:
            msg = 'Got wrong z index ({0}), not in [0,2]'
            print(msg.format(jl))
            return False
        if jx == jy:
            print('x and y indexes cannot be equal')
            return False
        if jx == jz:
            print('x and lambda indexes cannot be equal')
            return False
        if jy == jz:
            print('y and z indexes cannot be equal')
            return False

        # Define dimensions
        if self.__nx is None or self.__ny is None or \
           self.__nz is None:
            self.__nx = shape[jx]
            self.__ny = shape[jy]
            self.__nz = shape[jz]
        else:
            # Sanity check
            if not self.__sanity_check((shape[jx],shape[jy], \
                                        shape[jz])):
                return False

        # Save data
        self.__data[self.__idx_dic[var.lower()]] = \
                                            data.transpose((jx,jy,jz))

        # If input is logtau
        if var.lower() == 'ltau':
            self.__ltau = True
        elif var.lower() == 'tau':
            self.__ltau = False

        # Define precision
        if self.__prec is None:
            self.__prec = data.dtype

        # Success
        return True

######################################################################
######################################################################

    def create_file(self):
        ''' Create a file with the available data
        '''

        # A last sanity check
        if not self.__sanity_check():
            return False

        # Check necessary variables are here
        if not self.__check_vars():
            return False

        #
        # Try to create file
        try:

            # Open file
            f = open(self.__ofile,'wb')

            # Label
            f.write('2Dat'.encode())

            # Precision
            if self.__prec == np.float32:
                f.write(struct.pack('i',4))
            else:
                f.write(struct.pack('i',8))

            # Dimensions
            f.write(struct.pack('i',self.__nx))
            f.write(struct.pack('i',self.__ny))
            f.write(struct.pack('i',self.__nz))

            # Close
            f.close()

            # Let's save it via memmap
            data = np.memmap(self.__ofile, \
                             dtype=self.__prec, \
                             mode='r+', \
                             offset=20, \
                             shape=(self.__nx, \
                                    self.__ny, \
                                    24, \
                                    self.__nz))

            # Now fill variables
            for i in range(24):

                # If no variable, skip
                if self.__data[i] is None:

                    # Make zeros
                    data[:,:,i,:] = 0.
                    continue

                # If ltau
                if i==1 and self.__ltau:
                    # Asign
                    if self.__orderz:
                        data[:,:,i,:] = 10.**self.__data[i]
                    else:
                        data[:,:,i,:] = 10.**self.__data[i][:,:,::-1]

                # Rest
                else:

                    # Asign
                    if self.__orderz:
                        data[:,:,i,:] = self.__data[i]
                    else:
                        data[:,:,i,:] = self.__data[i][:,:,::-1]

            # Flush and free memmap
            data.flush()
            del data

        # Something failed
        except:
            raise

        # Success
        return True

######################################################################
######################################################################
######################################################################
######################################################################

class atmo1D_class():
    ''' Class to help with the creation of an atmospheric model for
        1D synthesis in HanleRT-TIC
    '''

######################################################################
######################################################################

    def __init__(self,afile,bfile=None):
        ''' Class initialization
        '''

        # Set output file path
        self.__afile = afile
        self.__bfile = bfile

        # Initialize variables to empty
        self.__data = []
        for i in range(16):
            self.__data.append(None)
        self.__flags = np.zeros((16),dtype=np.bool)
        self.__comm = None
        self.__name = None
        self.__wave = None

        self.__nz = None
        self.__h = None
        self.__ne = None
        self.__ltau = False

        # Index dictionary
        self.__idx_dic = {'h': 0, \
                          'tau': 0, \
                          'ltau': 0, \
                          't': 1, \
                          'pg': 2, \
                          'rho': 2, \
                          'bx': 13, \
                          'by': 14, \
                          'bz': 15, \
                          'vx': 5, \
                          'vy': 6, \
                          'vz': 3, \
                          'vmi': 4, \
                          'pe': 2, \
                          'rhoe': 2, \
                          'ne': 2, \
                          'nh0': 7, \
                          'nh1': 8, \
                          'nh2': 9, \
                          'nh3': 10, \
                          'nh4': 11, \
                          'np': 12}
        # Reverse index dictionary
        self.__idx_dic_r = {1: 't', \
                            13: 'bx', \
                            14: 'by', \
                            15: 'bz', \
                            5: 'vx', \
                            6: 'vy', \
                            3: 'vz', \
                            4: 'vmi', \
                            7: 'nh0', \
                            8: 'nh1', \
                            9: 'nh2', \
                            10: 'nh3', \
                            11: 'nh4', \
                            12: 'np'}

        # Descriptions dictionary
        self.__des_dic = {'h': 'height [km]', \
                          'tau/ltau': 'optical depth ' + \
                                      '(not/yes logarithmic)', \
                          't': 'temperature [K]', \
                          'pg': 'Gas pressure [dyn cm^-2]', \
                          'rho': 'Density [cm^-3]',  \
                          'bx': 'Magnetic field X component [G]', \
                          'by': 'Magnetic field Y component [G]', \
                          'bz': 'Magnetic field Z component [G]', \
                          'vx': 'Velocity X component [km s^-1]', \
                          'vy': 'Velocity Y component [km s^-1]', \
                          'vz': 'Velocity Z component [km s^-1]', \
                          'vmi': 'Microturbulent velocity', \
                          'pe': 'Electron pressure [dyn cm^-2]', \
                          'rhoe': 'Electron density [cm^-3]',  \
                          'ne': 'Electron number density [cm^-3]', \
                          'nht': 'Total hydrogen number ' + \
                                 'density [cm^-3]', \
                          'nha': 'Atomic hydrogen number ' + \
                                 'density [cm^-3]', \
                          'nhm': 'H- number density [cm^-3]', \
                          'nh0': 'Ground level hydrogen number ' + \
                                 'density [cm^-3]', \
                          'nh1': 'First excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'nh2': 'Second excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'nh3': 'Third excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'nh4': 'Fourth excited level hydrogen ' + \
                                 'number density [cm^-3]', \
                          'np': 'Proton number density [cm^-3]'}

######################################################################
######################################################################

    def __sanity_check(self,nz=None):
        ''' Check consistent dimensions
        '''

        # If specified nz
        if nz is not None:
            if self.__nz != nz:
                msg = 'The z dimensions ({0}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({1})'
                print(msg.format(nz,self.__nz))
                return False
            return True

        # For each variable
        for ivar,var in enumerate(self.__data):
            # Skip undefined
            if not self.__flags[ivar]: continue
            if ivar == 0:
                svar = self.__h
            elif ivar == 2:
                svar = self.__ne
            else:
                svar = self.__idx_dic_r[ivar]
            # Check dimensions
            if self.__nz != var.size:
                msg = 'The dimension for variable {0} ({1}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({2})'
                print(msg.format(svar,var.size,self.__nz))
                return False

        # Fine
        return True

######################################################################
######################################################################

    def __check_vars(self):
        ''' Check that all necessary variables are loaded
        '''
        if self.__data[0] is None:
            print('Either the height or the optical depth must be ' +\
                  'defined')
            return False
        else:
            diff = np.min(self.__data[0][1] - \
                          self.__data[0][0])
            if self.__h == 'h':
                if diff < 0.:
                    self.__order = True
                elif diff > 0:
                    self.__order = False
                else:
                    print('There cannot be two equal height values')
                    return False
            else:
                if diff > 0.:
                    self.__order = True
                elif diff < 0:
                    self.__order = False
                else:
                    print('There cannot be two equal optical depths')
                    return False
        if self.__data[1] is None:
            print('The temperature must be defined')
            return False
        if self.__data[2] is None:
            print('Either the gass pressure, the density, the ' + \
                  'electron pressure, the electron density, or ' + \
                  'the electron number density must be defined')
            return False
        return True

######################################################################
######################################################################

    def variable_list(self):
        ''' Return the list of available variables
        '''
        for key in self.__des_dic:
            print(f'- {key:3s}: {self.__des_dic[key]}')

######################################################################
######################################################################

    def set_comment(self,comm):
        ''' Define a comment
        '''
        # Sanity
        if not isinstance(comm,str):
            print('The comment must be a string')
            return False
        self.__comm = comm
        return True

######################################################################
######################################################################

    def set_name(self,name):
        ''' Define a name for the model
        '''
        # Sanity
        if not isinstance(name,str):
            print('The name must be a string')
            return False
        self.__name = name
        return True

######################################################################
######################################################################

    def set_wave(self,wave):
        ''' Define the wavelength
        '''
        # Sanity
        if not isinstance(comm,float) and \
           not isinstance(comm,int):
            print('The wavelength must be a number')
            return False
        self.__wave = wave
        return True

######################################################################
######################################################################

    def load_data(self,data,var):
        ''' Load the data in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The data must be a numpy array')
            return False
        # Sanity check
        if not isinstance(var, str):
            print('second argument must be a string')
            return False
        if var.lower() not in self.__idx_dic:
            print(f'Variable {var} not in the list, check with ' + \
                  'variable_list() method')
            return False

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims != 1:
            print('Unexpected number of dimensions in data. ' + \
                  'Got {0} instead of 1'.format(dims))
            return False

        # Define dimensions
        if self.__nz is None:
            self.__nz = shape[0]
        else:
            # Sanity check
            if not self.__sanity_check(shape[0]):
                return False

        # Save data
        self.__data[self.__idx_dic[var.lower()]] = data
        self.__flags[self.__idx_dic[var.lower()]] = True

        # If input is logtau
        if var.lower() == 'ltau':
            self.__ltau = True
        elif var.lower() == 'tau':
            self.__ltau = False

        # If index 0 or 2
        if self.__idx_dic[var.lower()] == 0:
            self.__h = var
        if self.__idx_dic[var.lower()] == 2:
            self.__ne = var

        # Success
        return True

######################################################################
######################################################################

    def create_file(self):
        ''' Create a file with the available data
        '''

        # A last sanity check
        if not self.__sanity_check():
            return False

        # Check necessary variables are here
        if not self.__check_vars():
            return False

        #
        # Try to create file
        try:

            # Open file for atmosphere
            f = open(self.__afile,'w')

            # Comment
            if self.__comm is None:

                # Generic
                f.write('* Model atmosphere created with ' + \
                        'standard tools\n')

            #
            else:

                # Length of comments
                comm = self.__comm
                words = comm.split()
                siz = len(' '.join(words))

                # While there are comments
                while siz > 0:

                    # Add first word
                    lcomm = words[0]
                    ii = 1

                    # While true
                    while True:
                        if ii >= len(words): break
                        if len(lcomm + ' ' + words[ii]) <= 80:
                            lcomm += ' ' + words[ii]
                            ii += 1
                        else:
                            break

                    # If we got something
                    f.write(f' * {lcomm}\n')

                    # Remove from words
                    if ii >= len(words):
                        siz = 0
                    else:
                        words = words[ii:]
                        siz = len(' '.join(words))

            # Name
            if self.__name is None:
                f.write('  MODEL\n')
            else:
                f.write(f'  {self.__name}\n')

            # Scale
            if self.__h == 'h':
                scal = '  HEIGHT SCALE'
                varz = 'HEIGHT [km]'
            else:
                scal = ' TAU SCALE'
                varz = '        TAU'
            if self.__wave is not None:
                scal += f' {self.__wave}'
            f.write(f'{scal}\n')

            # Gravity
            f.write('*\n* LG G\n  4.44\n')

            # Nodes
            f.write(f'*\n* NDEP\n {self.__nz}\n')

            # If horizontal velocities
            hvel = self.__flags[5] or self.__flags[6]

            # Header
            head = f'*\n*     {varz}'
            head += '         TEMP [K]'
            if self.__ne == 'pg':
                head += '  GAS PRES. [cgs]'
            if self.__ne == 'rho':
                head += '    DENSITY [cgs]'
            if self.__ne == 'ne':
                head += '       NE [cm^-3]'
            if self.__ne == 'pe':
                head += ' ELE. PRES. [cgs]'
            if self.__ne == 'rhoe':
                head += '  ELE DENS. [cgs]'
            head += '        VZ [km/s]'
            head += '   V_MICRO [km/s]'
            if hvel:
                head += '        VX [km/s]'
                head += '        VY [km/s]'
            head += '\n'
            f.write(head)

            # Get limits
            if self.__order:
                iz0 = 0
                iz1 = self.__nz
                diz = 1
            else:
                iz0 = self.__nz-1
                iz1 = -1
                diz = -1

            # Data
            for iz in range(iz0,iz1,diz):
                for ivar in range(5):
                    if ivar == 0:
                      if self.__h != 'h' and self.__ltau:
                        f.write(f'{10.**self.__data[ivar][iz]:17.8e}')
                      else:
                        f.write(f'{self.__data[ivar][iz]:17.8e}')
                    elif self.__flags[ivar]:
                        f.write(f'{self.__data[ivar][iz]:17.8e}')
                    else:
                        f.write(f'{0.:17.8e}')
                if hvel:
                    if self.__flags[5]:
                        f.write(f'{self.__data[5][iz]:17.8e}')
                    else:
                        f.write(f'{0.:17.8e}')
                    if self.__flags[6]:
                        f.write(f'{self.__data[6][iz]:17.8e}')
                    else:
                        f.write(f'{0.:17.8e}')
                f.write('\n')

            # If any other than 'ne'
            if self.__ne != 'ne':
                f.write('{self.__ne}\n')
            # If ne
            else:

                # Check if hydrogen
                hden = np.any(self.__flags[7:13])

                # If there is hydrogen
                if hden:

                    # Head
                    f.write('*\n* HYDROGEN POPULATIONS [cm^-3]\n')
                    f.write('*           NH(1)')
                    f.write('            NH(2)')
                    f.write('            NH(3)')
                    f.write('            NH(4)')
                    f.write('            NH(5)')
                    f.write('               NP\n')

                    # Heights
                    for iz in range(iz0,iz1,diz):
                        for i in range(7,13):
                            if self.__flags[i]:
                                f.write(f'{self.__data[i][iz]:17.8e}')
                            else:
                                f.write(f'{0.:17.8e}')
                        f.write('\n')

                # No hydrogen
                else:
                    f.write(f'{self.__ne}\n')
                

            # Close
            f.close()

            #
            # Magnetic field?
            if self.__bfile is not None:

                # Magnetic field?
                if np.any(self.__flags[13:16]):

                    # Open file for magnetic field
                    f = open(self.__bfile,'wb')

                    # Size and units
                    f.write(f'{self.__nz}\nDEG\n')

                    # Conversion
                    r2d = 180./np.pi

                    # For each height
                    for iz in range(iz0,iz1,diz):

                        # Bx
                        if self.__flags[13]:
                            Bx = self.__data[13][iz]
                        else:
                            Bx = 0e0
                        # By
                        if self.__flags[14]:
                            By = self.__data[14][iz]
                        else:
                            By = 0e0
                        # Bz
                        if self.__flags[15]:
                            Bz = self.__data[15][iz]
                        else:
                            Bz = 0e0

                        # Get polar
                        Bs = np.sqrt(Bx*Bx + By*By + Bz*Bz)

                        # If field
                        if Bs > 0.:

                            # Z component
                            if np.absolute(Bz) > 0.:

                                # Bh
                                Bh = np.sqrt(Bx*Bx + By*By)

                                # Horizontal components
                                if Bh > 0:

                                    # Inclination
                                    Bt = np.arccos(Bz/B)*r2d

                                    # No Y
                                    if np.absolute(By) <= 0.:
                                        if Bx > 0.:
                                            Bp = 0.
                                        else:
                                            Bp = 180.
                                    # No X
                                    elif np.absolute(Bx) <= 0.:
                                        if By > 0.:
                                            Bp = 90.
                                        else:
                                            Bp = 270.
                                    # Both X-Y
                                    else:
                                        Bp = np.arctan2(By,Bx)*r2d

                                # No horizontal
                                else:

                                    # Sign
                                    if Bz > 0.:
                                        Bt = 0.
                                    else:
                                        Bt = 180.
                                    Bp = 0.

                            # No z
                            else:

                                # Horizontal
                                Bt = 90.

                                # No Y
                                if np.absolute(By) <= 0.:
                                    if Bx > 0.:
                                        Bp = 0.
                                    else:
                                        Bp = 180.
                                # No X
                                elif np.absolute(Bx) <= 0.:
                                    if By > 0.:
                                        Bp = 90.
                                    else:
                                        Bp = 270.
                                # Both X-Y
                                else:
                                    Bp = np.arctan2(By,Bx)*r2d

                        # No field
                        else:
                            Bt = 0.
                            Bp = 0.

                        f.write(f'{Bs:16.8f}')
                        f.write(f'{Bt:16.8f}')
                        f.write(f'{Bp:16.8f}\n')

                    # Close file
                    f.close()

        # Something failed
        except:
            raise

        # Success
        return True

######################################################################
######################################################################
######################################################################
######################################################################

class weight_class():
    ''' Class to help with the creation of weight files for inversions
        in HanleRT-TIC
    '''

######################################################################
######################################################################

    def __init__(self,ofile):
        ''' Class initialization
        '''

        # Set output file path
        self.__ofile = ofile

        # Initialize variables to empty
        self.__data = []
        for i in range(4):
            self.__data.append(None)

        # Initialize default
        self.__default = 1e0

        # Size
        self.__nl = None

        # Index dictionary
        self.__idx_dic = {'I': 0, \
                          'Q': 1, \
                          'U': 2, \
                          'V': 3}

######################################################################
######################################################################

    def __sanity_check(self,nl=None):
        ''' Check consistent dimensions
        '''

        # If specified nz
        if nl is not None:
            if self.__nl != nl:
                msg = 'The size of the array ({0}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({1})'
                print(msg.format(nl,self.__nl))
                return False
            return True

        # For each variable
        for svar in list(self.__idx_dic):
            ivar = self.__idx_dic[svar]
            var = self.__data[ivar]
            # Skip undefined
            if var is None: continue
            # Check dimensions
            if self.__nl != var.size:
                msg = 'The dimension for variable {0} ({1}) ' + \
                      'is in conflict with the previously ' + \
                      'defined size ({2})'
                print(msg.format(svar,var.size,self.__nl))
                return False

        # Fine
        return True

######################################################################
######################################################################

    def set_default(self,val):
        ''' Set the default weight
        '''

        # Leave if no float
        if isinstance(val,np.ndarray):
            if val.size > 1:
                print('The input must be a number')
                return False
            else:
                self.__default = float(val[0])
                return True
        elif isinstance(val,float) or \
             isinstance(val,int):
            self.__default = float(val)
            return True
        else:
            print('The input must be a number')
            return False

######################################################################
######################################################################

    def variable_list(self):
        ''' Return the list of available variables
        '''
        for key in self.__idx_dic:
            print(f'- {key:3s}')

######################################################################
######################################################################

    def load_data(self,data,var):
        ''' Load the data in the class
        '''

        # Leave if no numpy
        if not isinstance(data, np.ndarray):
            print('The data must be a numpy array')
            return False
        # Sanity check
        if not isinstance(var, str):
            print('Second argument must be a string')
            return False
        # Sanity check
        if var not in list(self.__idx_dic):
            print(f'Variable {var} not found')
            return False

        # Get shape of the numpy array
        shape = data.shape
        dims = len(shape)

        # Unexpected shape
        if dims != 1:
            print('Unexpected number of dimensions in data. ' + \
                  'Got {0} instead of 1'.format(dims))
            return False

        # Define dimensions
        if self.__nl is None:
            self.__nl = shape[0]
        else:
            # Sanity check
            if not self.__sanity_check(shape[0]):
                return False

        # Save data
        self.__data[self.__idx_dic[var.upper()]] = data

        # Success
        return True

######################################################################
######################################################################

    def create_file(self):
        ''' Create a file with the available data
        '''

        # A last sanity check
        if not self.__sanity_check():
            return False

        #
        # Try to create file
        try:

            # Open file for atmosphere
            f = open(self.__ofile,'wb')

            # Write header
            f.write(struct.pack('i',self.__nl))
            f.write(struct.pack('i',4))

            # Close file
            f.close()

            # We write with memmap
            data = np.memmap(self.__ofile, \
                             dtype=np.float64, \
                             mode='r+', \
                             offset=8, \
                             shape=(4,self.__nl))

            # For each data
            for i,W in enumerate(self.__data):
                if W is None:
                    data[i,:] = self.__default
                else:
                    data[i,:] = W

            # Flush and free
            data.flush()
            del data

        # Something failed
        except:
            raise

        # Success
        return True

