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

