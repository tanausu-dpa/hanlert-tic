import os,struct
import numpy as np

def _verbose(msg):
    ''' Verbosity management
    '''
    # Verbose if allowed
    print(msg)

def _error(msg,level,ret=False):
    ''' Error management
    '''
    # Get prefix
    if level == 0:
        prefix = ' # Warning '
    elif level == 1:
        prefix = ' # Error '

    # Verbose
    _verbose(f'{prefix} {msg}')

    # If error and had to return something
    if ret and level == 1: return None

# Numpy integer
npint = type(np.argmin(np.array([0.])))

######################################################################
######################################################################
######################################################################

class _stokes_1D():
    ''' Class to manage emergent Stokes parameters from 1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_th': \
          [None,'Get LOS heliocentric angle'], \
         'get_ph': \
          [None,'Get LOS azimuthal angle'], \
         'get_mu': \
          [None,'Get LOS cosine of the heliocentric angle'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_stokesi': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get intensity [SI]'], \
          'get_stokesq': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get Stokes Q parameter'], \
          'get_stokesu': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get Stokes U parameter'], \
          'get_stokesv': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get Stokes V parameter'], \
          'get_linear': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get total linear polarization'], \
          'get_stokes': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get all Stokes parameters'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert emergence 1D file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__th = struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)
            self.__hsize = 22
            f.close()

            # Create memmaps
            self.__omg = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__hsize, \
                                   dtype=np.float64, \
                                   shape=(self.__nl))
            self.__stk = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__hsize + \
                                          8*self.__nl, \
                                   dtype=np.float64, \
                                   shape=(4,self.__nl))
            return True
        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def __get_gen_stokes(self,minl=None,maxl=None, \
                         fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters
        '''

        # Initialize
        out = [None,None,None,None]

        try:

            # Get wavelength
            if minl is not None or maxl is not None:
                lam = 1e2/self.__omg[::-1]

            # Intensity
            if 0 in indx or fractional:

                # Get intensity
                stkI = self.__stk[0,::-1]

                # Out?
                if 0 in indx:
                    out[0] = stkI.copy()

            # Q, U, and V
            for j in range(1,4):

                # If in output
                if j in indx:

                    # Read Stokes
                    out[j] = self.__stk[j,::-1].copy()

                # If fractional and output
                if fractional and j in indx:
                    out[j] /= stkI

            # Limits
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                for j in indx:
                    out[j] = out[j][i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                for j in indx:
                    out[j] = out[j][:i+1]

            # Norm
            for j in indx:
                if j == 0 or not fractional:
                    out[j] *= self.__unit_trans

            # Return
            return out

        except:
            raise

    def _get_stokesi(self,minl=None,maxl=None):
        ''' Get intensity from file
        '''
        return self.__get_gen_stokes(minl,maxl,False,[0])[0]

    def _get_stokesq(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes Q from file
        '''
        return self.__get_gen_stokes(minl,maxl,fractional,[1])[1]

    def _get_linear(self,minl=None,maxl=None,fractional=False):
        ''' Get total linear polarization from file
        '''
        qu = self.__get_gen_stokes(minl,maxl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokesu(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes U from file
        '''
        return self.__get_gen_stokes(minl,maxl,fractional,[2])[2]

    def _get_stokesv(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes V from file
        '''
        return self.__get_gen_stokes(minl,maxl,fractional,[3])[3]

    def _get_linear(self,minl=None,maxl=None,fractional=False):
        ''' Get total linear polarization from file
        '''
        qu = self.__get_gen_stokes(minl,maxl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes parameters from file
        '''
        iquv = self.__get_gen_stokes(minl,maxl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

######################################################################
######################################################################
######################################################################

class _stokesquad_1D():
    ''' Class to manage emergent Stokes parameters in the quadrature
        from 1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_nl': \
          [None,'Get number of directions'], \
         'get_th': \
          [None,'Get list of LOS heliocentric angle'], \
         'get_ph': \
          [None,'Get list of LOS azimuthal angle'], \
         'get_mu': \
          [None,'Get list of LOS cosine of the heliocentric angle'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_stokesi': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get intensity [SI]'], \
          'get_stokesq': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get Stokes Q parameter'], \
          'get_stokesu': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get Stokes U parameter'], \
          'get_stokesv': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get Stokes V parameter'], \
          'get_linear': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get total linear polarization'], \
          'get_stokes': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Get all Stokes parameters'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert emergence 1D file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__jump_to_lambda = 6
            f.seek(self.__nl*8,1)
            self.__nth = struct.unpack('i',f.read(4))[0]
            self.__nph = struct.unpack('i',f.read(4))[0]
            # Count dirs
            self.__nd = 0
            self.__th = []
            self.__ph = []
            for it in range(self.__nth):
              for ip in range(self.__nph):
                th = struct.unpack('d',f.read(8))[0]
                ph = struct.unpack('d',f.read(8))[0]
                f.seek(self.__nl*4*8,1)
                if th <= 90.0:
                    self.__nd += 1
                    self.__th.append(th)
                    self.__ph.append(ph)
            self.__th = np.array(self.__th)
            self.__ph = np.array(self.__ph)
            self.__mu = np.cos(self.__th*np.pi/180.)
            self.__hsize = 14 + 8*self.__nl
            f.close()

            # Create mamps
            self.__omg = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__jump_to_lambda, \
                                   dtype=np.float64, \
                                   shape=(self.__nl))
            self.__stk = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__hsize, \
                                   dtype=np.float64, \
                                   shape=(self.__nd,4,self.__nl))
            return True
        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_nd(self):
        ''' Get number of directions
        '''
        return self.__nd

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def __get_gen_stokes(self,minl=None,maxl=None, \
                         fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters
        '''

        # Initialize
        out = [None,None,None,None]

        # If cutting lambda
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        try:

            # Get space
            stk = np.empty((self.__nd,4,self.__nl))

            # For each direction
            for di in range(self.__nd):

                # Skip angles
                f.seek(16,1)

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = self.__stk[:,0,::-1]

                    # Out?
                    if 0 in indx:
                        out[0] = stkI.copy()

                # Q, U, and V
                for j in range(1,4):

                    # If in output
                    if j in indx:

                        # Read Stokes
                        out[j] = self.__stk[:,j,::-1].copy()

                    # If fractional and output
                    if fractional and j in indx:
                      out[j] /= out[j]/stkI

            # Limits
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                for j in indx:
                    out[j] = out[j][:,i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                for j in indx:
                    out[j] = out[j][:,:i+1]

            # Norm
            for j in indx:
                if j == 0 or not fractional:
                    out[j] *= self.__unit_trans

            # Return
            return out

        except:
            raise

    def _get_stokesi(self,minl=None,maxl=None):
        ''' Get intensity from file
        '''
        return self.__get_gen_stokes(minl,maxl,False,[0])[0]

    def _get_stokesq(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes Q from file
        '''
        return self.__get_gen_stokes(minl,maxl,fractional,[1])[1]

    def _get_linear(self,minl=None,maxl=None,fractional=False):
        ''' Get total linear polarization from file
        '''
        qu = self.__get_gen_stokes(minl,maxl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokesu(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes U from file
        '''
        return self.__get_gen_stokes(minl,maxl,fractional,[2])[2]

    def _get_stokesv(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes V from file
        '''
        return self.__get_gen_stokes(minl,maxl,fractional,[3])[3]

    def _get_linear(self,minl=None,maxl=None,fractional=False):
        ''' Get total linear polarization from file
        '''
        qu = self.__get_gen_stokes(minl,maxl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes(self,minl=None,maxl=None,fractional=False):
        ''' Get Stokes parameters from file
        '''
        iquv = self.__get_gen_stokes(minl,maxl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

######################################################################
######################################################################
######################################################################

class _contribution_1D():
    ''' Class to manage the contribution function from 1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e0

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_nz': \
          [None,'Get number of heights'], \
         'get_th': \
          [None,'Get LOS heliocentric angle'], \
         'get_ph': \
          [None,'Get LOS azimuthal angle'], \
         'get_mu': \
          [None,'Get LOS cosine of the heliocentric angle'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_height': \
          [{'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get vertical axis, either heights in [km] or ' + \
           'in optical depth'], \
          'get_ctrI': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get intensity contribution function [SI]'], \
          'get_ctrq': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get Stokes Q contribution function [SI]'], \
          'get_ctrU': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get Stokes U contribution function [SI]'], \
          'get_ctrV': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get Stokes V contribution function [SI]'], \
          'get_ctr': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get contribution function for all Stokes parameters'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D contribution function file head
        '''
        try:

            # Get head info
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__th = struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)
            self.__hsize = 26
            f.seek(self.__nl*8,1)
            z = struct.unpack('d'*self.__nz, \
                              f.read(8*self.__nz))
            self.__zreverse = z[-2] > z[-1]
            f.close()

            # Create memmaps
            self.__omg = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__hsize, \
                                   dtype=np.float64, \
                                   shape=(self.__nl))
            self.__z = np.memmap(self.__filename, \
                                 mode='r', \
                                 offset=self.__hsize + \
                                        self.__nl*8, \
                                 dtype=np.float64, \
                                 shape=(self.__nz))
            self.__ctr = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__hsize + \
                                          self.__nl*8 +  \
                                          self.__nz*8, \
                                   dtype=np.float64, \
                                   shape=(4,self.__nz,self.__nl))
            return True

        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def _get_height(self,minh=None,maxh=None):
        ''' Get height from file
        '''
        if self.__zreverse:
            iminh = maxh
            imaxh = minh
        else:
            iminh = minh
            imaxh = maxh
        try:
            z = self.__z.copy()
            if self.__zreverse:
                z = np.log10(z)
            else:
                z *= 1e-5
            if iminh is not None:
                i = np.argmin(np.absolute(z - iminh))
                z = z[i:]
            if imaxh is not None:
                i = np.argmin(np.absolute(z - imaxh))
                z = z[:i+1]
            return z
        except:
            raise

    def __get_gen_ctr(self,minl=None,maxl=None, \
                      minh=None,maxh=None,indx=[0]):
        ''' Generic read of contribution function
        '''

        # Reversed z?
        if self.__zreverse:
            iminh = maxh
            imaxh = minh
        else:
            iminh = minh
            imaxh = maxh

        # Initialize
        out = [None,None,None,None]

        # Get data
        try:

            # Cutting lambda
            if minl is not None or maxl is not None:
                lam = 1e2/self.__omg[::-1]

            # Cutting height
            if minh is not None or maxh is not None:
                z = self.__z.copy()
                if self.__zreverse:
                    z = np.log10(z)
                else:
                    z *= 1e-5

            # For each Stokes
            for j in range(4):

                # Read
                if j in indx:

                    out[j] = self.__ctr[j,:,::-1].copy()

            # Cut
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                for j in indx:
                    out[j] = out[j][:,i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                for j in indx:
                    out[j] = out[j][:,:i+1]
            if minh is not None:
                i = np.argmin(np.absolute(z - minh))
                z = z[i:]
                for j in indx:
                    out[j] = out[j][i:,:]
            if maxh is not None:
                i = np.argmin(np.absolute(z - maxh))
                z = z[:i+1]
                for j in indx:
                    out[j] = out[j][:i+1,:]

            # Units
            for j in indx:
                out[j] *= self.__unit_trans

            # Return
            return out

        except:
            raise

    def _get_ctri(self,minl=None,maxl=None,minh=None,maxh=None):
        ''' Get intensity contribution function from file
        '''
        return self.__get_gen_ctr(minl,maxl,minh,maxh,[0])[0]

    def _get_ctrq(self,minl=None,maxl=None,minh=None,maxh=None):
        ''' Get Stokes Q contribution function from file
        '''
        return self.__get_gen_ctr(minl,maxl,minh,maxh,[1])[1]

    def _get_ctru(self,minl=None,maxl=None,minh=None,maxh=None):
        ''' Get Stokes U contribution function from file
        '''
        return self.__get_gen_ctr(minl,maxl,minh,maxh,[2])[2]

    def _get_ctrv(self,minl=None,maxl=None,minh=None,maxh=None):
        ''' Get Stokes V contribution function from file
        '''
        return self.__get_gen_ctr(minl,maxl,minh,maxh,[3])[3]

    def _get_ctr(self,minl=None,maxl=None,minh=None,maxh=None):
        ''' Get all contribution function from file
        '''
        iquv =  self.__get_gen_ctr(minl,maxl,minh,maxh,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

######################################################################
######################################################################
######################################################################

class _tau_1D():
    ''' Class to manage the height for optical depth equal to 1 from
        1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_th': \
          [None,'Get LOS heliocentric angle'], \
         'get_ph': \
          [None,'Get LOS azimuthal angle'], \
         'get_mu': \
          [None,'Get LOS cosine of the heliocentric angle'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_height': \
          [{'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get height for optical depth equal to 1 in [km] or ' + \
           'in optical depth'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D height where optical depth is one file
            head
        '''
        try:

            # Get header
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__th = struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)
            self.__hsize = 26
            f.close()

            # Create memmaps
            self.__omg = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__hsize, \
                                   dtype=np.float64, \
                                   shape=(self.__nl))
            self.__tau1 = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__hsize + \
                                           8*self.__nl, \
                                    dtype=np.float64, \
                                    shape=(self.__nl))

            return True
        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def _get_height(self,minl=None,maxl=None):
        ''' Get height for optical depth equal to 1
        '''
        try:

            if minl is not None or maxl is not None:
                lam = 1e2/self.__omg[::-1]

            tau1 = self.__tau1[::-1].copy()*1e-5
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                tau1 = tau1[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                tau1 = tau1[:i+1]
            return tau1
        except:
            raise

######################################################################
######################################################################
######################################################################

class _jkq_1D():
    ''' Class to manage the frequency integrated radiation field
        tensors from a 1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_na': \
          [None,'Get number of atoms'], \
         'get_nz': \
          [None,'Get number of height'], \
         'get_nt': \
          [{'atom': 'Particular atom index to get number of ' + \
                    'transitions'}, \
           'Get number of transitions'], \
          'get_height': \
          [{'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get vertical axis, either heights in [km] or \
            in optical depth'], \
          'get_jkq': \
          [{'atom': 'Request this particular atom index', \
            'transition': 'Request this particular transition ' + \
                          'index. Requires "atom"', \
            'k': 'Request this particular multipole. ' + \
                 'Requires "transition"', \
            'q': 'Request this particular multipolar component. ' + \
                 'Requires "k"', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height', \
            'sti': \
             'Request radiation field tensor integrated over the ' + \
             'emissivity, if available'}, \
           'Get the frequency integrated radiation field tensors'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D JKQ file head
        '''
        try:

            # Get actual header
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__sti = bool(struct.unpack('?', f.read(1))[0])
            f.seek(3,1) # bool padding
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__na = struct.unpack('i',f.read(4))[0]
            self.__nxtran = struct.unpack('i',f.read(4))[0]
            self.__hsize = 18

            # Check z
            z = struct.unpack('d'*self.__nz, \
                              f.read(8*self.__nz))
            self.__zreverse = z[-2] > z[-1]

            # Save atomic data and jump size
            self.__nt = np.zeros((self.__na),dtype=int)
            self.__jump = np.zeros((self.__na),dtype=int)

            # For each atom
            for ia in range(self.__na):
                # Read transitions
                self.__nt[ia] = struct.unpack('i',f.read(4))[0]
                self.__jump[ia] = 18*8*self.__nt[ia]*self.__nz
                f.seek(self.__jump[ia],1)

            # Close file
            f.close()

            # Create memmaps
            self.__z = np.memmap(self.__filename, \
                                 mode='r', \
                                 offset=self.__hsize, \
                                 dtype=np.float64, \
                                 shape=(self.__nz))

            return True
        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_na(self):
        ''' Get number of atoms
        '''
        return self.__na

    def _get_nt(self,atom=None):
        ''' Get number of transitions
        '''
        try:
            return self.__nt[atom]
        except:
            return self.__nt

    def _get_height(self,minh=None,maxh=None):
        ''' Get height from file
        '''
        if self.__zreverse:
            iminh = maxh
            imaxh = minh
        else:
            iminh = minh
            imaxh = maxh
        try:
            z = self.__z.copy()
            if self.__zreverse:
                z = np.log10(z)
            else:
                z *= 1e-5
            if iminh is not None:
                i = np.argmin(np.absolute(z - iminh))
                z = z[i:]
            if imaxh is not None:
                i = np.argmin(np.absolute(z - imaxh))
                z = z[:i+1]
            return z
        except:
            raise

    def _get_jkq(self,atom=None,transition=None,k=None,q=None, \
                 minh=None,maxh=None,sti=False):
        ''' Get frequency integrated radiation field tensors
        '''

        # Manage height limits
        if self.__zreverse:
            iminh = maxh
            imaxh = minh
        else:
            iminh = minh
            imaxh = maxh

        # Check it is possible
        if sti and not self.__sti:
            _error('No stimulated emission in the file',0,True)
            return None

        # Try fetching data
        try:

            # If cutting
            if minh is not None or maxh is not None:
                z = self.__z.copy()
                if self.__zreverse:
                    z = np.log10(z)
                else:
                    z *= 1e-5

            # Open file
            f = open(self.__filename,'rb')
            f.seek(self.__hsize + self.__nz*8,0)

            # No atom specified
            if atom is None:

                # Big output
                jkq = np.zeros((np.sum(self.__nt),3,5,2,self.__nz))

                # Running index
                ind = -1

                # For each atom
                for ia in range(self.__na):

                    # If asking for stimulated, skip normal
                    if sti:
                        f.seek(np.sum(self.__jump[ia]+4),1)
                    # Skip transition
                    else:
                        f.seek(4,1)

                    # For each transition
                    for it in range(self.__nt[ia]):

                        # Run index
                        ind += 1

                        # For each K
                        for K in range(3):
                            for iQ,Q in enumerate(range(-K,K+1)):

                                # Adjust iQ
                                jQ = iQ + 2 - K

                                # Read
                                jkq[ind,K,jQ,:,:] = np.transpose( \
                             np.array( \
                              struct.unpack('d'*2*self.__nz, \
                                            f.read(8*2*self.__nz))). \
                              reshape((self.__nz,2)), (1,0))

                # Adjust height
                if iminh is not None:
                    i = np.argmin(np.absolute(z - iminh))
                    z = z[i:]
                    jkq = jkq[:,:,:,:,i:]
                if imaxh is not None:
                    i = np.argmin(np.absolute(z - imaxh))
                    z = z[:i+1]
                    jkq = jkq[:,:,:,:,:i+1]


            # Atom specified
            else:

                # Jump to atom
                if atom > 0:

                    # Jump to current atom
                    f.seek(np.sum(self.__jump[:atom]) + 4*atom,1)

                    # Skip normal if stimulated
                    if sti:
                        f.seek(np.sum(self.__jump[atom]),1)

                # Transition not specified
                if transition is None:

                    # Output
                    jkq = np.zeros((self.__nt[atom],3,5,2,self.__nz))

                    # For each transition
                    for it in range(self.__nt[atom]):

                        # For each K
                        for K in range(3):
                            for iQ,Q in enumerate(range(-K,K+1)):

                                # Adjust iQ
                                jQ = iQ + 2 - K

                                # Read
                                jkq[it,K,jQ,:,:] = np.transpose( \
                            np.array( \
                              struct.unpack('d'*2*self.__nz, \
                                            f.read(8*2*self.__nz))). \
                              reshape((self.__nz,2)), (1,0))

                    # Adjust height
                    if iminh is not None:
                        i = np.argmin(np.absolute(z - iminh))
                        z = z[i:]
                        jkq = jkq[:,:,:,:,i:]
                    if imaxh is not None:
                        i = np.argmin(np.absolute(z - imaxh))
                        z = z[:i+1]
                        jkq = jkq[:,:,:,:,:i+1]

                # Transition  specified
                else:

                    # Output
                    jkq = np.zeros((3,5,2,self.__nz))

                    # Jump to transition
                    f.seek(transition*18*8*self.__nz,1)

                    # For each K
                    for K in range(3):
                        for iQ,Q in enumerate(range(-K,K+1)):

                            # Adjust iQ
                            jQ = iQ + 2 - K

                            # Read
                            jkq[K,jQ,:,:] = np.transpose( \
                             np.array( \
                              struct.unpack('d'*2*self.__nz, \
                                            f.read(8*2*self.__nz))). \
                              reshape((self.__nz,2)), (1,0))

                    # Adjust height
                    if iminh is not None:
                        i = np.argmin(np.absolute(z - iminh))
                        z = z[i:]
                        jkq = jkq[:,:,:,i:]
                    if imaxh is not None:
                        i = np.argmin(np.absolute(z - imaxh))
                        z = z[:i+1]
                        jkq = jkq[:,:,:,:i+1]

                    # If particular K
                    if k is not None:

                        # Trim
                        jkq = jkq[k,:,:,:]

                        # If particular Q
                        if q is not None:

                            # Get index
                            jq = q + 2

                            # Trim
                            jkq = jkq[jq,:,:]
            f.close()
            return jkq*self.__unit_trans
        except struct.error:
            raise
        except:
            raise

######################################################################
######################################################################
######################################################################

class _rkq_1D():
    ''' Class to manage the density matrix from a 1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_na': \
          [None,'Get number of atoms'], \
         'get_nz': \
          [None,'Get number of height'], \
         'get_height': \
          [{'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get vertical axis, either heights in [km] ', + \
           'or in optical depth'], \
          'get_rkq': \
          [{'atom': 'Request this particular atom index', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height'}, \
           'Get the density matrix tensor components'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D JKQ file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__na = struct.unpack('i',f.read(4))[0]
            self.__hsize = 10

            # Check z
            z = np.array(struct.unpack('d'*self.__nz, \
                                       f.read(8*self.__nz)))
            self.__zreverse = z[-2] > z[-1]

            # Close file
            f.close()

            # Create memmaps
            self.__z = np.memmap(self.__filename, \
                                 mode='r', \
                                 offset=self.__hsize, \
                                 dtype=np.float64, \
                                 shape=(self.__nz))

            return True
        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_na(self):
        ''' Get number of atoms
        '''
        return self.__na

    def _get_height(self,minh=None,maxh=None):
        ''' Get height from file
        '''
        if self.__zreverse:
            iminh = maxh
            imaxh = minh
        else:
            iminh = minh
            imaxh = maxh
        try:
            z = self.__z.copy()
            if self.__zreverse:
                z = np.log10(z)
            else:
                z *= 1e-5
            if iminh is not None:
                i = np.argmin(np.absolute(z - iminh))
                z = z[i:]
            if imaxh is not None:
                i = np.argmin(np.absolute(z - imaxh))
                z = z[:i+1]
            return z
        except:
            raise

    def _get_rkq(self,atom=None,minh=None,maxh=None):
        ''' Get frequency integrated radiation field tensors
        '''

        # Manage height limits
        if self.__zreverse:
            iminh = maxh
            imaxh = minh
        else:
            iminh = minh
            imaxh = maxh

         # Manage atom input
        if atom is None:
            iatoms = list(range(self.__na))
        else:
            if not isinstance(atom,list):
                _error('atom argument must be a list',0,True)
                return None
            for ia in atom:
                if not isinstance(ia,int):
                    _error('atom argument must be a list ' + \
                           'of integers',0,True)
                    return None
            iatoms = copy.deepcopy(atom)

        # Try fetching data
        try:

            # If cutting
            if minh is not None or maxh is not None:
                z = self.__z.copy()
                if self.__zreverse:
                    z = np.log10(z)
                else:
                    z *= 1e-5

            # Open file
            f = open(self.__filename,'rb')
            f.seek(self.__hsize+8*self.__nz,0)

            # Adjust height
            if iminh is not None:
                i0 = np.argmin(np.absolute(z - iminh))
            else:
                i0 = 0
            if imaxh is not None:
                i1 = np.argmin(np.absolute(z - imaxh))
            else:
                i1 = self.__nz

            # Initialize
            rkq = []

            # For each atom
            for ia in range(self.__na):

                # Check if storing
                keep = ia in iatoms

                # Initialize
                if keep: lrkq = {}

                # Population
                if keep:
                    lrkq['n'] = \
                          np.array(struct.unpack('d'*self.__nz, \
                                                 f.read(8*self.__nz)))
                else:
                    f.seek(self.__nz*8,1)

                # Number of terms
                nt = struct.unpack('i',f.read(4))[0]
                if keep: lrkq['nt'] = nt

                # For each term
                for it in range(nt):

                    # Initialize
                    if keep: lrkq[it] = {}

                    # J levels
                    nj = struct.unpack('i',f.read(4))[0]
                    if keep: lrkq[it]['nj'] = nj
                    njj = nj*nj
                    if keep: lrkq[it]['J'] = []

                    # Multi-term
                    mt =  nj > 1

                    # For each J combination
                    for jj in range(njj):

                        # Read J combination
                        J1 = struct.unpack('i',f.read(4))[0]
                        if keep:
                            if J1 not in lrkq[it]['J']:
                                lrkq[it]['J'].append(J1)
                        J2 = struct.unpack('i',f.read(4))[0]
                        if keep:
                            if J2 not in lrkq[it]['J']:
                                lrkq[it]['J'].append(J2)

                        # Get label
                        if keep:
                            if J1 % 2 == 0:
                                ijs1 = f'{int(round(J1*0.5))}'
                            else:
                                ijs1 = f'{J1}/2'
                            if ijs1 not in lrkq[it]:
                                lrkq[it][ijs1] = {}
                            if mt:
                                if J2 % 2 == 0:
                                    ijs2 = f'{int(round(J2*0.5))}'
                                else:
                                    ijs2 = f'{J2}/2'
                                lrkq[it][ijs1][ijs2] = {}

                        # K limits
                        kmin = int(np.absolute(int(round(J1-J2))//2))
                        kmax = int(np.absolute(int(round(J1+J2))//2))

                        # Pointer
                        if keep:
                            if mt:
                                point = lrkq[it][ijs1][ijs2]
                            else:
                                point = lrkq[it][ijs1]

                        # For each K
                        for K in range(kmin,kmax+1):

                            # Initialize
                            if keep:
                                point[K] = {}

                            # For each Q
                            for Q in range(-K,K+1):

                                # Read
                                if keep:
                                    point[K][Q] = \
                                           np.zeros((self.__nz), \
                                                    dtype=np.complex_)
                                    for iz in range(self.__nz):
                                        point[K][Q][iz] = \
                                   struct.unpack('d',f.read(8))[0] + \
                                   1j*struct.unpack('d',f.read(8))[0]
                                        f.seek(4,1)
                                    point[K][Q] = point[K][Q][i0:i1+1]

                                # Skip
                                else:
                                    f.seek(self.__nz*20,1)

                    # Process J
                    if keep:
                        for jj in range(len(lrkq[it]['J'])):
                            lrkq[it]['J'][jj] = 0.5*lrkq[it]['J'][jj]
                        lrkq[it]['J'] = np.array(lrkq[it]['J'])

                # Add to list
                if keep:
                    rkq.append(lrkq.copy())
                    del lrkq

            f.close()
            return rkq
        except struct.error:
            raise
        except:
            raise

######################################################################
######################################################################
######################################################################

class _jkqnu_1D():
    ''' Class to manage the frequency dependent radiation field
        tensors from a 1D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        self.__methods = {\
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nz': \
          [None,'Get number of height'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_jkq': \
          [{'minl': \
             'Lower boundary for output height', \
            'maxl': \
             'Upper boundary for output height'}, \
           'Get the frequency dependent radiation field tensors'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D JKQnu file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__nl = struct.unpack('i',f.read(4))[0]
            f.close()

            # Sizes info
            self.__jump_to_lambda = 10
            self.__hsize = self.__jump_to_lambda + 8*self.__nl

            # File size
            real_size = os.path.getsize(self.__filename)

            # Check if anisotropy
            self.__anis =  real_size >= self.__hsize + \
                                        self.__nz*self.__nl*18*8

            # Create memmaps
            self.__omg = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__jump_to_lambda, \
                                   dtype=np.float64, \
                                   shape=(self.__nl))
            if self.__anis:
                self.__jkq = np.memmap(self.__filename, \
                                       mode='r', \
                                       offset=self.__hsize, \
                                       dtype=np.float64, \
                                       shape=(self.__nz, \
                                              self.__nl, \
                                              9,2))[:,::-1,:,:]
            else:
                self.__jkq = np.memmap(self.__filename, \
                                       mode='r', \
                                       offset=self.__hsize, \
                                       dtype=np.float64, \
                                       shape=(self.__nz, \
                                              self.__nl))[:,::-1]
            return True
        except struct.error:
            raise
        except:
            raise

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def _get_jkq(self,minl=None,maxl=None):
        ''' Get frequency integrated radiation field tensors
        '''

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        try:

            # Fetch
            data = self.__jkq.copy()

            # Transform
            data *= self.__unit_trans

            # Polarized
            if self.__anis:

                # Adjust wavelength
                if minl is not None:
                    i = np.argmin(np.absolute(lam - minl))
                    lam = lam[i:]
                    data = data[:,i:,:,:]
                if maxl is not None:
                    i = np.argmin(np.absolute(lam - maxl))
                    lam = lam[:i+1]
                    data = data[:,:i+1,:,:]

                # Get output
                jkq = {0: {0: data[:,:,0,0]}, \
                       1: {-1: data[:,:,1,0] + 1j*data[:,:,1,1], \
                           0: data[:,:,2,0], \
                           1: data[:,:,3,0] + 1j*data[:,:,3,1]}, \
                       2: {-2: data[:,:,4,0] + 1j*data[:,:,4,1], \
                           -1: data[:,:,5,0] + 1j*data[:,:,5,1], \
                           0: data[:,:,6,0], \
                           1: data[:,:,7,0] + 1j*data[:,:,7,1], \
                           2: data[:,:,8,0] + 1j*data[:,:,8,1]}}
                # Return
                return jkq

            # Unpolarized
            else:

                # Adjust wavelength
                if minl is not None:
                    i = np.argmin(np.absolute(lam - minl))
                    lam = lam[i:]
                    data = data[:,i:]
                if maxl is not None:
                    i = np.argmin(np.absolute(lam - maxl))
                    lam = lam[:i+1]
                    data = data[:,:i+1]

                # Get sizes
                nl = data.shape[1]
                nz = self.__nz

                # Get output
                jkq = {0: {0: data}, \
                       1: {-1: np.zeros((nz,nl)), \
                           0: np.zeros((nz,nl)), \
                           1: np.zeros((nz,nl))}, \
                       2: {-2: np.zeros((nz,nl)), \
                           -1: np.zeros((nz,nl)), \
                           0: np.zeros((nz,nl)), \
                           1: np.zeros((nz,nl)), \
                           2: np.zeros((nz,nl))}}
                # Return
                return jkq
        except:
            raise

######################################################################
######################################################################
######################################################################

class _cols_damp_1D():
    ''' Class to manage the 1D collisions, damping parameter, and
        elastic rates files
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
       #self.__unit_trans = 1e6
        #  Remain cgs
        self.__unit_trans = 1e0

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere']}
        if self.__cols:
            self.__methods['get_type'] = \
              [None,'Get the type of collisional rates in the file']
            self.__methods['get_nl'] = \
              [None,'Get number of levels or terms']
            self.__methods['get_dims'] = \
              [None,' Get number of positions in height axis, ' + \
                    'atoms, and levels/terms']
            self.__methods['get_data'] = \
            [{'ia': 'List of atoms to include in the output', \
              'i1': 'Initial level/term in the output', \
              'i2': 'Final level/term in the output'}, \
            'Extract the collisional rates']
        else:
            self.__methods['get_nt'] = \
              [None,'Get number of transitions']
            self.__methods['get_dims'] = \
              [None,' Get number of positions in height axis, ' + \
                    'atoms, and transitions']
            if self.__damp:
                self.__methods['get_data'] = \
                  [{'ia': 'List of atoms to include in the output', \
                    'it': 'Transitions in the output'}, \
                'Extract the damping parameter']
            else:
                self.__methods['get_data'] = \
                  [{'ia': 'List of atoms to include in the output', \
                    'it': 'Transition in the output'}, \
                'Extract the elastic rates']

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert population/departure coeff. file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')

            # Skip 1 first letter
            f.seek(1,0)

            # Read the second for the type of file
            self.__type = f.read(1).decode('utf-8')
            self.__cols = self.__type == 'l' or self.__type == 't'
            if not self.__cols:
                self.__damp = self.__type == 'a'

            # Get dimensions
            self.__na = int(struct.unpack('i',f.read(4))[0])
            self.__nz = int(struct.unpack('i',f.read(4))[0])

            # Initialize entries
            self.__nn = []
            self.__siz = []

            # For each atom
            for ia in range(self.__na):

                # Append size
                self.__nn.append(int(struct.unpack('i',f.read(4))[0]))

                # Collisions
                if self.__cols:
                    self.__siz.append(int(8.*self.__nz* \
                                          self.__nn[-1]* \
                                          self.__nn[-1]))
                else:
                    self.__siz.append(int(8.*self.__nz*self.__nn[-1]))

            # Make arrays
            self.__nn = np.array(self.__nn)
            self.__siz = np.array(self.__siz)

            # Close file
            f.close()

            # Size of head
            self.__head = 10

        except struct.error:
            return False
            raise
        except:
            raise

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_type(self):
        ''' Get type of collisions in this file
        '''
        if self.__type == 'l':
            _error('Collisions between level',0)
        elif self.__type == 't':
            _error('Collisions between terms',0)
        else:
            _error('Type not recognized',1)

    def _get_nl(self):
        ''' Get number of levels in file
        '''
        return self.__nn

    def _get_nt(self,atom=None):
        ''' Get number of transitions
        '''
        try:
            return self.__nn[atom]
        except:
            return self.__nn

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_dims(self):
        ''' Get number of positions in height axes, atoms, and levels,
            terms, or transitions
        '''
        return self.__nz, self.__na, self.__nn

    def _get_data(self,ia=None,i1=None,i2=None,it=None):
        ''' Get collisional rates/dampings for a given pixel
        '''

        # If no particular atom requested
        if ia is None:

            # No trimming data
            trim = False

            # Potential errors
            if self.__cols:
                if i1 is not None or i2 is not None:
                    _error('Must specify ia to specify i1 or i2',1)
            else:
                if it is not None:
                    _error('Must specify ia to specify it',1)

        # If selecting an atom
        else:

            # Trimming data
            trim = True

            # Generate a list if it is not
            if not isinstance(ia, list):
                ivar = [ia]
            else:
                ivar = ia.copy()

            # Check integers and bounded
            for avar in ivar:
                if not isinstance(avar, int):
                    _error('The field ia requires integers',1)
                    return None
                if avar < 0 or avar >= self.__na:
                    _error('The requested atoms ' + avar + \
                           ' are out of limits, ' + \
                           'check with get_nlevel',1)
                    return None

            # Error handling with collisional requests
            if len(ivar) > 1:
                if self.__cols:
                    if i1 is not None or i2 is not None:
                        _error('Must specify only one ia ' + \
                               'to specify i1 or i2',1)
                else:
                    if it is not None:
                        _error('Must specify only one ia ' + \
                               'to specify it',1)

        # If the file has collisions
        if self.__cols:

            # Check if requested origin
            if i1 is not None:
                if not isinstance(i1, int):
                    _error('The field i1 requires an integer',1)
                    return None
                if i1 < 0 or i1 >= self.__nn:
                    _error('The requested level ' + i1 + \
                           ' is out of limits, ' + \
                           'check with get_nl',1)
                    return None

            # Check if requested destiniy
            if i2 is not None:
                if not isinstance(i2, int):
                    _error('The field i2 requires an integer',1)
                    return None
                if i2 < 0 or i2 >= self.__nn:
                    _error('The requested level ' + i2 + \
                           ' is out of limits, ' + \
                           'check with get_nl',1)
                    return None

        # If damping parameters
        else:

            # Check if requested transition
            if it is not None:
                if not isinstance(it, int):
                    _error('The field it requires an integer',1)
                    return None
                if it < 0 or it >= self.__nn:
                    _error('The requested transition ' + it + \
                           ' is out of limits, ' + \
                           'check with get_nt',1)
                    return None

        # Outputs
        out = {}

        # Try getting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__head,0)

            # For each atom
            for ia in range(self.__na):

                # Skip
                f.seek(4,1)

                # Reading atom?
                if trim:
                    if ia not in ivar:
                        f.seek(self.__siz[ia],1)
                        continue

                # Read data
                col = struct.unpack('d'*(self.__siz[ia]//8), \
                                    f.read(self.__siz[ia]))

                # Process
                if self.__cols:
                    out[ia] = np.array(col).reshape((self.__nn[ia], \
                                                     self.__nn[ia], \
                                                     self.__nz))*1e8
                else:
                    out[ia] = np.array(col).reshape((self.__nn[ia], \
                                                     self.__nz))
            # Close
            f.close()

        except:
            raise

        # Check if trimming levels
        if trim:
            # Only one atom allowed
            if len(ivar) == 1:
                # If collisions
                if self.__cols:
                    # Based on arguments
                    if i1 is not None or i2 is not None:
                        if i1 is not None and i2 is not None:
                            out = out[ivar[0]][i1,i2,:]
                        elif i1 is not None:
                            out = out[ivar[0]][i1,:,:]
                        elif i2 is not None:
                            out = out[ivar[0]][:,i2,:]
                # Damping
                else:
                    if it is not None:
                        out = out[ivar[0]][it,:,:]

        # Return column
        return out

######################################################################
######################################################################
######################################################################

class _back_1D():
    ''' Class to manage the 1D background files
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e3
        #  Transformation to CGS
       #self.__unit_trans = 1e0/299792458e2

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nd': \
          [None,'Get number of directions in the model atmosphere'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the height, direction, ' + \
                'and wavelength dimensions'], \
         'get_vars': \
          [None,'Get list of continuum variables'], \
         'get_vars_alias': \
          [None,'Get list of continuum variables with their ' + \
                'aliases'], \
         'get_vars_units': \
          [None,'Get list of continuum variables with their ' + \
                'corresponding units'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_data': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the ' + \
             'available ones with get_vars_alias()}'}, \
           'Extract the continuum variables at a particular column']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert background continuum file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            f.seek(self.__nl*8,1)
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__nd = struct.unpack('i',f.read(4))[0]
            f.close()

        except struct.error:
            return False
            raise
        except:
            raise

        # Variables
        self.__nvar = 3
        self.__vars = [r'$\eta_{\rm c}$',r'$\sigma_{\rm c}$', \
                       r'$\epsilon{\rm c}$']
        self.__alias = ['eta','sig','eps']
        self.__vars_units = ['[m$^{-1}$]','[m$^{-1}$]', \
                             '[J m$^{-3}$ s$^{-1}$' + \
                             ' sr$^{-1}$ Hz$^{-1}$]']

        # Sizes
        self.__jump_to_lambda = 6
        self.__head = self.__jump_to_lambda + self.__nl*8 + 8
        self.__siz = self.__nz*self.__nl*self.__nd*self.__nvar

        # Create memmaps
        self.__omg = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        self.__data = np.memmap(self.__filename, \
                                mode='r', \
                                offset=self.__head, \
                                dtype=np.float64, \
                                shape=(self.__nz, \
                                       self.__nd, \
                                       self.__nvar,\
                                       self.__nl))

        # Return valid
        return True


    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nd(self):
        ''' Get number of directions
        '''
        return self.__nd

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in height, directions, and
            wavelength axes
        '''
        return self.__nz, self.__nd, self.__nl

    def _get_vars(self):
        ''' Get variables with node data
        '''
        return self.__vars

    def _get_vars_alias(self):
        ''' Get variables and their alias
        '''
        out = []
        for var,alias in zip(self.__vars,self.__alias):
            out.append(var+' -> ',alias)
        return out

    def _get_vars_units(self):
        ''' Get variables with node data with units
        '''
        out = []
        for var,uni in zip(self.__vars,self.__vars_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def _get_data(self,minl=None,maxl=None,var=None):
        ''' Get background quantities for a given pixel
        '''

        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # Try getting data
        try:

            # Fetch data
            col = self.__data[:,:,:,::-1].copy()

            # Adjust wavelength
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                col = col[:,:,:,i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                col = col[:,:,:,:i+1]

        except:
            raise

        # Return column
        out = {}
        if 'eta' in ivar:
            out['eta'] = col[:,:,0,:]
        if 'sig' in ivar:
            out['sig'] = col[:,:,1,:]
        if 'eps' in ivar:
            out['eps'] = col[:,:,2,:]*self.__unit_trans

        return out

######################################################################
######################################################################
######################################################################

class _pop_dep_1D():
    ''' Class to manage the 1D population/departure files
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
       #self.__unit_trans = 1e6
        #  Remain cgs
        self.__unit_trans = 1e0

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nl': \
          [None,'Get number of levels'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_dims': \
          [None,' Get number of positions in height axis ' + \
                'and levels'], \
         'get_data': \
          [{'ie': \
            'List of levels to include in the output'}, \
            'Extract the populations/departures']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert population/departure coeff. file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')

            # Skip 1 first letter
            f.seek(1,0)

            # Read the second for the type of file
            self.__type = f.read(1).decode('utf-8')

            # Get dimensions
            self.__nz = int(struct.unpack('i',f.read(4))[0])
            self.__nn = int(struct.unpack('i',f.read(4))[0])

            # Close file
            f.close()

            # Size of head
            self.__head = 10

            # Create memmap
            self.__pop = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__head, \
                                   dtype=np.float64, \
                                   shape=(self.__nz, \
                                          self.__nn))

        except struct.error:
            return False
            raise
        except:
            raise

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_type(self):
        ''' Get type of population/departure file
        '''
        if self.__type == 'p':
            _error('Number density',0)
        elif self.__type == 'b':
            _error('Departure coefficient',0)
        else:
            _error('Type not recognized',1)

    def _get_nl(self):
        ''' Get number of levels in file
        '''
        return self.__nn

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_dims(self):
        ''' Get number of positions in x, y, and height axes
            and entries
        '''
        return self.__nz, self.__nn

    def _get_data(self,ie=None):
        ''' Get populations/departures for a given pixel
        '''

        # If var is not None
        if ie is None:
            trim = False
        else:
            trim = True
            if not isinstance(ie, list):
                ivar = [ie]
            else:
                ivar = ie.copy()
            for evar in ivar:
                if not isinstance(evar, int):
                    _error('The field var requires integers',1)
                    return None
                if evar < 0 or evar >= self.__nn:
                    _error('The requested levels ' + evar + \
                           ' is out of limits, ' + \
                           'check with get_nlevel',1)
                    return None

        # Try getting data
        try:

            # Read data column
            col = self.__pop.copy()

        except:
            raise

        # If trimming
        if trim:
            todel = []
            for i in range(self.__nn-1,-1,-1):
                if i not in ivar:
                    todel.append(i)
            col = np.delete(col,np.array(todel,dtype='int32'),axis=1)

        # Units
        if self.__type == 'p':
            col *= self.__unit_trans

        # Return column
        return col

######################################################################
######################################################################
######################################################################

class _stokes_15D():
    ''' Class to manage emergent Stokes parameters from 15D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        # Method
        self.__methods = {
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_polarized': \
          [None,'Get if the synthesis included polarization'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, and ' + \
                'wavelength dimensions'], \
         'get_th': \
          [None,'Get LOS heliocentric angle'], \
         'get_ph': \
          [None,'Get LOS azimuthal angle'], \
         'get_mu': \
          [None,'Get LOS cosine of the heliocentric angle'], \
         'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
         'get_stokesi_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract Stokes I at a particular column'], \
         'get_stokesq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular column'], \
         'get_stokesu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular column'], \
         'get_stokesv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular column'], \
         'get_linear_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a ' + \
           'particular column'], \
         'get_stokes_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular column'], \
          'get_plane_stk': \
          [{'il': \
             'Coordinate in the wavelength dimension of the ' + \
             'Stokes parameters to extract', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract Stokes parameters at a given wavelength ' + \
           'position for ' + \
           'the whole field of view'], \
         'get_stokesi_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract Stokes I at a particular wavelength ' + \
           'index for ' + \
           'the whole field of view'], \
         'get_stokesq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_linear_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a ' + \
           'particular wavelength index for the whole ' + \
           'field of view'], \
         'get_stokes_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular ' + \
           'wavelength index for the whole field of view'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method. Note that ' + \
           'the raw data is not scaled to the SI units and ' + \
          f'you need to multiply by {self.__unit_trans}. ' + \
           'You can get this number with the method get_scale()'], \
          'get_scale': \
          [None,f'Get the scale factor between the raw data and ' + \
           'SI units'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D emergence file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)
            # Wavelength
            self.__nl = struct.unpack('i',f.read(4))[0]
            # Dimensions
            self.__nx = struct.unpack('i',f.read(4))[0]
            self.__ny = struct.unpack('i',f.read(4))[0]
            # LOS
            self.__th= struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)
            # Lambda
            self.__jump_to_lambda = 4*4 + 8*2
            # Head
            self.__head = self.__jump_to_lambda + self.__nl*8

        except struct.error:
            raise
        except:
            raise

        # Close file
        f.close()

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size
        expectedIsize = self.__head + self.__nl*self.__nx* \
                                      self.__ny*4
        expectedPsize = self.__head + self.__nl*self.__nx* \
                                      self.__ny*4*4

        # If intensity size
        if real_size == expectedIsize:

            self.__mode = 1
            self.__complete = True
            self.__column_size = self.__nl*4

        # If polarization size
        elif real_size == expectedPsize:

            self.__mode = 2
            self.__complete = True
            self.__column_size = self.__nl*16

        # Incomplete file
        else:

            # Try reading four columns worth of Stokes intensity
            try:

                # Open
                f = open(self.__filename,'rb')

                # Put in head position
                f.seek(self.__head,0)

                # Read what could be polarization column
                test = np.array(struct.unpack('f'*self.__nl*4, \
                                              f.read(4*self.__nl*4)))

                # And close
                f.close()

                # Test sign
                if real_size > expectedIsize or \
                   np.min(test) < 0.:

                    # Guessed polarization
                    self.__mode = 2
                    self.__complete = False
                    msg = f'I have guessed that this is ' + \
                          f'an incomplete polarization file'
                    _error(msg,0)
                    self.__column_size = self.__nl*16

                # No sign change
                else:

                    # Guessed intensity
                    self.__mode = 1
                    self.__complete = False
                    msg = f'I have guessed that this is ' + \
                          f'an incomplete intensity file'
                    _error(msg,0)
                    self.__column_size = self.__nl*4

            except:

                # Failed
                msg = f'Expected size {expectedIsize} or ' + \
                      f'{expectedPsize}, but got {real_size} ' + \
                      f'instead'
                _error(msg,0)
                raise

        # Create memmaps
        self.__omg = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        # If complete
        if self.__complete:
            # Intensity
            if self.__mode == 1:
                self.__data = np.memmap(self.__filename, \
                                        mode='r', \
                                        offset=self.__head, \
                                        dtype=np.float32, \
                                        shape=(self.__nx, \
                                               self.__ny, \
                                               1,self.__nl))
            # Full Stokes
            else:
                self.__data = np.memmap(self.__filename, \
                                        mode='r', \
                                        offset=self.__head, \
                                        dtype=np.float32, \
                                        shape=(self.__nx, \
                                               self.__ny, \
                                               4,self.__nl))

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_polarized(self):
        ''' Get if there is polarization
        '''
        return self.__mode > 1

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in x, y, and wavelength axes
        '''
        return self.__nx, self.__ny, self.__nl

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_scale(self):
        ''' Return the multiplicative factor to get SI units
        '''
        return self.__unit_trans

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def __get_gen_column(self,ix,iy,minl=None,maxl=None, \
                         fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters column
        '''

        # Initialize output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # If the file is complete
        if self.__complete:

            # Try getting data
            try:

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = self.__data[ix,iy,0,::-1].copy()

                    # Out?
                    if 0 in indx:
                        out[0] = stkI

                # Q, U, and V
                for j in range(1,4):

                    # There is polarization
                    if self.__mode == 2:

                        # To output
                        if j in indx:

                            # Read Stokes
                            out[j] = self.__data[ix,iy,j,::-1].copy()

                    # No pol
                    else:

                        # Zeros
                        out[j] = np.zeros((self.__nl))

                    # Manage units
                    if fractional and j in indx:
                        out[j] /= stkI

            # Others
            except:
                raise

        # If the file is incomplete
        else:

            # Get size to read
            siz = self.__nl
            bsiz = siz*4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*self.__column_size + \
                           self.__ny*ix*self.__column_size,0)

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = np.array(struct.unpack('f'*siz, \
                                                  f.read(bsiz)))[::-1]

                    # Out?
                    if 0 in indx:
                        out[0] = stkI

                # No intensity
                else:

                    # Skip
                    f.seek(bsiz,1)

                # Q, U, and V
                for j in range(1,4):

                    # There is polarization
                    if self.__mode == 2:

                        # To output
                        if j in indx:

                            # Read Stokes
                            out[j] = np.array( \
                                    struct.unpack('f'*siz, \
                                                  f.read(bsiz)))[::-1]
                        else:

                            # Skip
                            f.seek(bsiz,1)

                    # No pol
                    else:

                        # Zeros
                        out[j] = np.zeros((siz))

                    # Manage units
                    if fractional and j in indx:
                        out[j] /= stkI

            # Failed
            except struct.error:

                # If the file is complete, the error is more severe,
                # let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Warn
                    msg = 'Could not read, may be due to ' + \
                          'the file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    for j in indx:
                        out[j] = np.zeros((self.__nl))

            # Others
            except:
                raise

            # Close file
            f.close()

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            for j in indx:
                out[j] = out[j][i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            for j in indx:
                out[j] = out[j][:i+1]

        # Units
        if fractional:
            if 0 in indx:
                out[0] *= self.__unit_trans
        else:
            for j in indx:
                out[j] *= self.__unit_trans

        # Return
        return out

    def _get_stokesi_column(self,ix,iy,minl=None,maxl=None):
        ''' Get intensity profile at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl,False,[0])[0]

    def _get_stokesq_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes Q profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[1])[1]

    def _get_stokesu_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes U profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[2])[2]

    def _get_stokesv_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes V profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix,int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy,int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[3])[3]

    def _get_linear_column(self,ix,iy,minl=None,maxl=None, \
                           fractional=False):
        ''' Get Stokes V profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix,int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy,int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        qu = self.__get_gen_column(ix,iy,minl,maxl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_column(self,ix,iy,minl=None,maxl=None, \
                           fractional=False):
        ''' Get Stokes parameter at a given column
        '''

        # Mode?
        if self.__mode == 1:
           _error('The file is only intensity',1)
           return None

        # Valid?
        if not isinstance(ix,int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy,int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        iquv = self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def __get_gen_plane(self,il,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters plane
        '''

        # Output
        out = [None,None,None,None]

        # If file is complete
        if self.__complete:

            try:

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = self.__data[:,:,0,il].copy()

                    # Out?
                    if 0 in indx:
                        out[0] = stkI

                # Q, U, and V
                for j in range(1,4):

                    # There is polarization
                    if self.__mode == 2:

                        # To output
                        if j in indx:

                            # Get Stokes
                            out[j] = self.__data[:,:,j,il].copy()

                    # No polarization
                    else:

                        # Zero
                        out[j] = np.zeros((self.__nx,self.__ny))

                    # Manage units
                    if fractional and j in indx:
                        out[j][ix,iy] /= stkI

            except:
                raise

        # Incomplete file
        else:

            # Get size to read
            left = il*4
            right = (self.__nl - il - 1)*4
            full = self.__nl*4
            abort = False

            # Output
            out = [None,None,None,None]

            # For each index requested
            for j in indx:
                out[j] = np.empty((self.__nx,self.__ny))

            # Open file
            f = open(self.__filename,'rb')

            # Seek to data
            f.seek(self.__head,0)

            # For each column
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Intensity
                        if 0 in indx or fractional:

                            # Get intensity
                            if left > 0: f.seek(left,1)
                            stkI = struct.unpack('f',f.read(4))[0]
                            if right > 0: f.seek(right,1)

                            # Out?
                            if 0 in indx:
                                out[0][ix,iy] = stkI

                        # No intensity
                        else:

                            # Skip
                            f.seek(full,1)

                        # Q, U, and V
                        for j in range(1,4):

                            # There is polarization
                            if self.__mode == 2:

                                # To output
                                if j in indx:

                                    # Get Stokes
                                    if left > 0: f.seek(left,1)
                                    out[j][ix,iy] = \
                                       struct.unpack('f',f.read(4))[0]
                                    if right > 0: f.seek(right,1)

                                # No output
                                else:

                                    # Skip
                                    f.seek(full,1)

                            # No polarization
                            else:

                                # Zero
                                out[j][ix,iy] = 0.0

                            # Manage units
                            if fractional and j in indx:
                                out[j][ix,iy] /= stkI

                    # Reading error
                    except struct.error:

                        # If the file is complete, the error is
                        # more severe, let it crash
                        if self.__complete:
                            raise

                        # Incomplete file, may be missing data
                        else:

                            # Warn
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)

                            # Generate zeros
                            for j in indx:
                                out[j][ix,iy:self.__ny] = 0.0
                            abort = True
                            break

                    except:
                        raise

                # We are leaving
                if abort:
                    for j in indx:
                        out[j][ix+1:self.__nx,:] = 0.0
                    break

            # Close file
            f.close()

        # Units
        if fractional:
            if 0 in indx:
                out[0] *= self.__unit_trans
        else:
            for j in indx:
                out[j] *= self.__unit_trans

        # Return
        return out

    def _get_stokesi_plane(self,il):
        ''' Get intensity profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,False,[0])[0]

    def _get_stokesq_plane(self,il,fractional=False):
        ''' Get Stokes Q profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,fractional,[1])[1]

    def _get_stokesu_plane(self,il,fractional=False):
        ''' Get Stokes U profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,fractional,[2])[2]

    def _get_stokesv_plane(self,il,fractional=False):
        ''' Get Stokes V profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,fractional,[3])[3]

    def _get_linear_plane(self,il,fractional=False):
        ''' Get Stokes V profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        qu = self.__get_gen_plane(jl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_plane(self,il,fractional=False):
        ''' Get Stokes profiles at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        iquv = self.__get_gen_plane(jl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def _get_cube(self):
        ''' Get Stokes profiles memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__data[:,:,:,::-1]

######################################################################
######################################################################
######################################################################

class _contribution_15D():
    ''' Class to manage the contribution function from 15D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e0

        # Method
        self.__methods = {
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height ' + \
                'dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, height, and ' + \
                'wavelength dimensions'], \
         'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_ctri_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get intensity contribution function [SI] at a ' + \
           'given column'], \
          'get_ctrq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to ' + \
             'extract', \
            'iy': \
             'Coordinate in the y dimension of the column to ' + \
             'extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get Stokes Q contribution function [SI] at a ' + \
           'given column'], \
          'get_ctru_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get Stokes U contribution function [SI] at a ' + \
           'given column'], \
          'get_ctrv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get Stokes V contribution function [SI] at a ' + \
           'given column'], \
          'get_ctr_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get full Stokes vector contribution function [SI]' + \
           ' at a given column'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method. Note that ' + \
           'the raw data is not scaled to the SI units and ' + \
          f'you need to multiply by {self.__unit_trans}. ' + \
           'You can get this number with the method get_scale()'], \
          'get_scale': \
          [None,f'Get the scale factor between the raw data and ' + \
           'SI units'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D contribution function file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Wavelength
            self.__nl = struct.unpack('i',f.read(4))[0]

            # Dimensions
            self.__nx = struct.unpack('i',f.read(4))[0]
            self.__ny = struct.unpack('i',f.read(4))[0]
            self.__nz = struct.unpack('i',f.read(4))[0]

            # LOS
            self.__th= struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)

            # Lambda
            self.__jump_to_lambda = 4*5 + 8*2

            # Head
            self.__head = self.__jump_to_lambda + self.__nl*8

        except struct.error:
            raise
        except:
            raise

        # Close file
        f.close()

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size
        expectedsize = self.__head + self.__nl*self.__nx* \
                                     self.__ny*self.__nz*4*4

        # Column size
        self.__column_size = self.__nl*self.__nz*4*4

        # If complete
        if real_size == expectedsize:

            self.__complete = True

        # Incomplete file
        else:

            # Incomplete file
            self.__complete = False
            msg = f'This is an incomplete file. \n' + \
                  f'Expected size {expectedsize}, ' + \
                  f'but got {real_size} instead'
            _error(msg,0)

        # Create memmaps
        self.__omg = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        if self.__complete:
            self.__data = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__head, \
                                    dtype=np.float32, \
                                    shape=(self.__nx, \
                                           self.__ny, \
                                           4,self.__nz, \
                                           self.__nl))

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nxyz(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self.__nx, self.__ny, self.__nz

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in x, y, height, and
            wavelength axes
        '''
        return self.__nx, self.__ny, self.__nz, self.__nl

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_scale(self):
        ''' Return the multiplicative factor to get SI units
        '''
        return self.__unit_trans

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def __get_gen_column(self,ix,iy,minl=None,maxl=None,indx=[0]):
        ''' Generic read of contribution function column
        '''

        # Initialize output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # If complete file
        if self.__complete:

            # Try getting data
            try:

                # For each stokes parameter
                for j in range(4):

                    # To output
                    if j in indx:

                        # Read Stokes
                        out[j] = self.__data[ix,iy,j,:,::-1].copy()

            # Others
            except:
                raise

        # Incomplete file
        else:

            # Get size to read
            siz = self.__nl*self.__nz
            bsiz = siz*4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*self.__column_size + \
                           self.__ny*ix*self.__column_size,0)

                # For each stokes parameter
                for j in range(4):

                    # To output
                    if j in indx:

                        # Read Stokes
                        out[j] = (np.array( \
                                  struct.unpack('f'*siz, \
                                                f.read(bsiz))). \
                                   reshape((self.__nz, \
                                            self.__nl)))[:,::-1]
                    else:

                        # Skip
                        f.seek(bsiz,1)

            # Failed
            except struct.error:

                # If the file is complete, the error is more severe,
                # let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Warn
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    for j in indx:
                        out[j] = np.zeros((self.__nz,self.__nl))

            # Others
            except:
                raise

            # Close file
            f.close()

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            for j in indx:
                out[j] = out[j][:,i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            for j in indx:
                out[j] = out[j][:,:i+1]

        # Units
        for j in indx:
            out[j] *= self.__unit_trans

        # Return
        return out

    def _get_ctri_column(self,ix,iy,minl=None,maxl=None):
        ''' Get intensity contribution function at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl,[0])[0]

    def _get_ctrq_column(self,ix,iy,minl=None,maxl=None):
        ''' Get Stokes Q contribution function at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl,[1])[1]

    def _get_ctru_column(self,ix,iy,minl=None,maxl=None):
        ''' Get Stokes U contribution function at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl,[2])[2]

    def _get_ctrv_column(self,ix,iy,minl=None,maxl=None):
        ''' Get Stokes V contribution function at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl,[3])[3]

    def _get_ctr_column(self,ix,iy,minl=None,maxl=None):
        ''' Get full Stokes vector contribution function at a
            given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        iquv =  self.__get_gen_column(ix,iy,minl,maxl,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def _get_cube(self):
        ''' Get contribution memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__data[:,:,:,:,::-1]

######################################################################
######################################################################
######################################################################

class _tau_15D():
    ''' Class to manage the height for optical depth equal to 1 from
        1.5D synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = {
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, and ' + \
                'wavelength dimensions'], \
         'get_th': \
          [None,'Get LOS heliocentric angle'], \
         'get_ph': \
          [None,'Get LOS azimuthal angle'], \
         'get_mu': \
          [None,'Get LOS cosine of the heliocentric angle'], \
         'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
         'get_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get height for optical depth equal to 1 in [km] or ' + \
           'in optical depth for a given column'], \
         'get_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Get height for optical depth equal to 1 in [km] or ' + \
           'in optical depth at a particular wavelength ' + \
           'index for the whole field of view'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method.'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D height where optical depth is one
            file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Wavelength
            self.__nl = struct.unpack('i',f.read(4))[0]

            # Dimensions
            self.__nx = struct.unpack('i',f.read(4))[0]
            self.__ny = struct.unpack('i',f.read(4))[0]

            # LOS
            self.__th= struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)

            # Lambda
            self.__jump_to_lambda = 4*4 + 8*2

            # Head
            self.__head = self.__jump_to_lambda + self.__nl*8

        except struct.error:
            raise
        except:
            raise

        # Close file
        f.close()

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size
        expectedsize = self.__head + self.__nl*self.__nx* \
                                     self.__ny*4

        # If intensity size
        if real_size == expectedsize:

            self.__complete = True

        # Incomplete file
        else:

            # Incomplete file
            self.__complete = False
            msg = f'This is an incomplete file. \n' + \
                  f'Expected size {expectedsize}, ' + \
                  f'but got {real_size} instead'
            _error(msg,0)

        # Create memmaps
        self.__omg = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        if self.__complete:
            self.__data = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__head, \
                                    dtype=np.float32, \
                                    shape=(self.__nx, \
                                           self.__ny, \
                                           self.__nl))

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in x, y, and wavelength axes
        '''
        return self.__nx, self.__ny, self.__nl

    def _get_th(self):
        ''' Get LOS heliocentric angle
        '''
        return self.__th

    def _get_ph(self):
        ''' Get LOS azimuthal angle
        '''
        return self.__ph

    def _get_mu(self):
        ''' Get LOS cosine of the heliocentric angle
        '''
        return self.__mu

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def _get_column(self,ix,iy,minl=None,maxl=None):
        ''' Get height for optical depth equal to 1 at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # If complete file
        if self.__complete:

            # Try getting data
            try:

                # Get intensity
                tau = self.__data[ix,iy,::-1].copy()

            # Others
            except:
                raise

        # Incomplete file
        else:

            # Get size to read
            siz = self.__nl
            bsiz = siz*4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*bsiz + self.__ny*ix*bsiz,0)

                # Get intensity
                tau = np.array(struct.unpack('f'*siz, \
                                              f.read(bsiz)))[::-1]

            # Failed
            except struct.error:

                # If the file is complete, the error is more severe,
                # let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Warn
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    tau = np.zeros((self.__nl))

            # Others
            except:
                raise

            # Close file
            f.close()

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            tau = tau[i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            tau = tau[:i+1]

        # Return
        return tau

    def _get_plane(self,jl):
        ''' Get height for optical depth equal to 1 for a given
            wavelength index
        '''

        # Valid?
        if not isinstance(jl, int) and not isinstance(jl, npint):
           _error('il must be an integer',1)
           return None
        if jl < 0 or jl >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        # Reverse index
        il = self.__nl - jl - 1

        # If complete file
        if self.__complete:

            # Output
            tau = self.__data[:,:,il]

        # Incomplete file
        else:

            # Get size to read
            left = il*4
            right = (self.__nl - il - 1)*4
            full = self.__nl*4
            abort = False

            # Output
            tau = np.empty((self.__nx,self.__ny))

            # Open file
            f = open(self.__filename,'rb')

            # Seek to data
            f.seek(self.__head,0)

            # For each column
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Get Stokes
                        if left > 0: f.seek(left,1)
                        tau[ix,iy] = struct.unpack('f',f.read(4))[0]
                        if right > 0: f.seek(right,1)

                    # Reading error
                    except struct.error:

                        # If the file is complete, the error
                        # is more severe, let it crash
                        if self.__complete:
                            raise

                        # Incomplete file, may be missing data
                        else:

                            # Warn
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)

                            # Generate zeros
                            tau[ix,iy:self.__ny] = 0.0
                            abort = True
                            break

                    except:
                        raise

                # We are leaving
                if abort:
                    tau[ix+1:self.__nx,:] = 0.0
                    break

            # Close file
            f.close()

        # Return
        return tau

    def _get_cube(self):
        ''' Get Stokes profiles memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__data[:,:,::-1]

######################################################################
######################################################################
######################################################################

class _atmo_1D():
    ''' Class to manage the 1D model atmosphere format
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_vars': \
          [None,'Get list of variables in the model atmosphere'], \
         'get_vars_units': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding units'], \
         'get_vars_alias': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding alias'], \
          'get_name': \
          [None,'Get the name of the model specified in the file'], \
          'get_comment': \
          [None,'Get the comment specified in the file'], \
          'get_wavelength': \
          [None,'Get the wavelength corresponding to the optical ' + \
           'depth scale present or to be calculated'], \
          'get_column': \
          [{'minz': \
             'Lower height or optical depth boundary (lower ' + \
             'limit) for output', \
            'maxz': \
             'Upper height or optical depth boundary (upper ' + \
             'limit) for output', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_alias()}'}, \
           'Extract the model atmosphere']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D model atmosphere
        '''
        try:

            # Initialize comment and variables
            self.__comment = ''
            self.__vars = []
            self.__alias = []
            self.__vars_units = []
            self.__i0 = 0
            self.__i1 = 0
            self.__block1 = []
            self.__block2 = []

            # Number of valid lines counter
            ival = 0

            # Open file
            f = open(self.__filename,'r')

            # For each line
            for iline,line in enumerate(f):

                # Strip comments
                if '!' in line:
                    il = line.find('!')
                    if ival < 1: cc = line[il+1:]
                    cline = line[:il]
                elif '*' in line:
                    il = line.find('*')
                    if ival < 1: cc = line[il+1:]
                    cline = line[:il]
                else:
                    cline = line
                    cc = ''

                # Strip spaces
                cline = cline.strip()

                # Comment
                if ival < 1:
                    cc = cc.strip()
                    if len(cc) > 1:
                        if self.__comment == '':
                            self.__comment = cc
                        else:
                            self.__comment += ' '+cc

                # If empty, skip
                if len(cline) < 1: continue

                # There is content
                ival += 1

                # The first one is the name
                if ival == 1: self.__name = cline

                # The second valid line is the scale
                # or tau
                if ival == 2:
                    if 'height' in cline.lower():
                        self.__ztype = 0
                    elif 'tau' in cline.lower():
                        self.__ztype = 1
                    try:
                        self.__wavelength = float(cline.split()[-1])
                    except:
                        self.__wavelength = 500.

                # The third # valid line is a float
                if ival == 3: self.__logg = float(cline)

                # The fourth valid line is the size
                if ival == 4: self.__nz = int(cline)

                # The fifth is the first line of the block, we
                # can check if we have horizontal velocities
                if ival == 5:

                    # Save line
                    self.__i0 = iline

                    # Check existence
                    self.__vx = len(cline.split()) > 5
                    self.__vy = len(cline.split()) > 6

                # Skip if too early
                if ival < 5: continue

                # The nz + 5 can be either a label or a list
                # of numbers, telling us the type
                if ival == 5+self.__nz:

                    # Save
                    self.__i1 = iline

                    # 6 Columns?
                    if len(cline.split()) == 6:
                        nh = True
                        for c in cline.split():
                            try:
                                c = float(c)
                            except:
                                nh = False
                                break
                    else:
                        nh = False
                    if nh:
                        self.__etype = 0
                    elif cline.lower() == 'ne':
                        self.__etype = 1
                    elif cline.lower() == 'pg':
                        self.__etype = 2
                    elif cline.lower() == 'rho':
                        self.__etype = 3
                    elif cline.lower() == 'pe':
                        self.__etype = 4
                    elif cline.lower() == 'rhoe':
                        self.__etype = 5
                    else:
                        self.__etype = -1


            #
            # We can now fill the variables completely

            # Vertical axis
            if self.__ztype == 0:
                self.__vars.append('h')
                self.__alias.append('h')
                self.__vars_units.append('km')
            else:
                self.__vars.append(r'$\log{\tau_{\rm c}}$')
                self.__alias.append('ltau')
                self.__vars_units.append('')

            # Temperature
            self.__vars.append('T')
            self.__alias.append('T')
            self.__vars_units.append('K')

            # NE variable
            if self.__etype == 0 or self.__etype == 1:
                self.__vars.append(r'N$_{\rm e}$')
                self.__alias.append('ne')
                self.__vars_units.append(r'cm$^{-3}$')
            elif self.__etype == 2:
                self.__vars.append(r'P$_{\rm g}$')
                self.__alias.append('Pg')
                self.__vars_units.append(r'dyn cm$^{-2}$')
            elif self.__etype == 3:
                self.__vars.append(r'$\rho$')
                self.__alias.append('rho')
                self.__vars_units.append(r'g cm$^{-3}$')
            elif self.__etype == 4:
                self.__vars.append(r'P$_{\rm e}$')
                self.__alias.append('Pe')
                self.__vars_units.append(r'dyn cm$^{-2}$')
            elif self.__etype == 5:
                self.__vars.append(r'$\rho_{\rm e}$')
                self.__alias.append('rhoe')
                self.__vars_units.append(r'g cm$^{-3}$')

            # Velocity
            if self.__vx:
                self.__vars.append(r'v$_{\rm x}')
                self.__alias.append('vx')
                self.__vars_units.append(r'km s$^{-1}$')
            if self.__vy:
                self.__vars.append(r'v$_{\rm y}')
                self.__alias.append('vy')
                self.__vars_units.append(r'km s$^{-1}$')
            self.__vars.append(r'v$_{\rm z}')
            self.__alias.append('vz')
            self.__vars_units.append(r'km s$^{-1}$')
            self.__vars.append(r'v$_{\rm mi}')
            self.__alias.append('vmi')
            self.__vars_units.append(r'km s$^{-1}$')

            # These variables are block 1
            for var in self.__alias:
                self.__block1.append(var)

            # Hydrogen
            if self.__etype == 0:
                self.__vars.append(r'N$_{\rm H_0}$')
                self.__vars.append(r'N$_{\rm H_1}$')
                self.__vars.append(r'N$_{\rm H_2}$')
                self.__vars.append(r'N$_{\rm H_3}$')
                self.__vars.append(r'N$_{\rm H_4}$')
                self.__vars.append(r'N$_{\rm p^+}$')
                self__alias.append('nH0')
                self__alias.append('nH1')
                self__alias.append('nH2')
                self__alias.append('nH3')
                self__alias.append('nH4')
                self__alias.append('np')
                self.__vars_units.append(r'cm$^{-3}$')
                self.__vars_units.append(r'cm$^{-3}$')
                self.__vars_units.append(r'cm$^{-3}$')
                self.__vars_units.append(r'cm$^{-3}$')
                self.__vars_units.append(r'cm$^{-3}$')
                self.__vars_units.append(r'cm$^{-3}$')
                self__block2.append('nH0')
                self__block2.append('nH1')
                self__block2.append('nH2')
                self__block2.append('nH3')
                self__block2.append('nH4')
                self__block2.append('np')

        except:
            raise

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_name(self):
        ''' Get the name in the read file
        '''
        return self.__name

    def _get_comment(self):
        ''' Get the comment in the read file
        '''
        return self.__name

    def _get_wavelength(self):
        ''' Get wavelength for the continuum corresponding to
            the optical depth
        '''
        return self.__wavelength

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_vars(self):
        ''' Get variables in model atmosphere
        '''
        return self.__vars

    def _get_vars_alias(self):
        ''' Get variables and their alias in model atmosphere
        '''
        out = []
        for var,alias in zip(self.__vars,self.__alias):
            out.append(var+' -> '+alias)
        return out

    def _get_vars_units(self):
        ''' Get variables in model atmosphere with units
        '''
        out = []
        for var,uni in zip(self.__vars,self.__vars_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_column(self,minz=None,maxz=None,var=None):
        ''' Get model atmosphere
        '''

        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Initialize variables
        out = {}
        for var in ivar:
            out[var] = np.zeros((self.__nz))

        # Asking for block1?
        block1 = False
        for var in ivar:
            if var in self.__block1:
                block1 = True
                break

        # Asking for block2?
        block2 = False
        for var in ivar:
            if var in self.__block2:
                block2 = True
                break

        # Initialize
        needz = 0
        needne = False


        # If asked for ne variable
        for var in ivar:
            if 'h' == var:
                needz = 1
                zvar = var
            elif 'ltau' == var:
                needz = 1
                zvar = var
            if 'ne' == var:
                needne =True
                nevar = var
            elif 'pe' == var:
                needne = True
                nevar = var
            elif 'rhoe' == var:
                needne = True
                nevar = var
            elif 'pg' == var:
                needne = True
                nevar = var
            elif 'rho' == var:
                needne = True
                nevar = var

        # Complete
        if minz is not None or maxz is not None: needz += 10


        # Try getting data
        try:

            # Open file
            with open(self.__filename,'r') as f:

                # Run over lines
                for iline, line in enumerate(f):

                    # If need block 1
                    if block1:

                        # If in the block
                        if iline >= self.__i0 and \
                           iline < self.__i0+self.__nz:

                            # Strip
                            if '!' in line:
                                idx = line.find('!')
                                cline = line[:idx].strip()
                            elif '*' in line:
                                idx = line.find('*')
                                cline = line[:idx].strip()
                            else:
                                cline = line.strip()

                            # Cols
                            cols = cline.split()

                            # Height index
                            iz = iline - self.__i0

                            if needz > 0:
                                out[zvar][iz] = float(cols[0])
                            if 'T' in ivar:
                                out['T'][iz] = float(cols[1])
                            if needne:
                                out[nevar][iz] = float(cols[2])
                            if 'vz' in ivar:
                                out['vz'][iz] = float(cols[3])
                            if 'vmi' in ivar:
                                out['vmi'][iz] = float(cols[4])
                            if 'vx' in ivar:
                                out['vx'][iz] = float(cols[5])
                            if 'vy' in ivar:
                                out['vy'][iz] = float(cols[6])

                    # If need block 2
                    if block2:

                        # If in the block
                        if iline >= self.__i1 and \
                           iline < self.__i1+self.__nz:

                            # Height index
                            iz = iline - self.__i0

                            # Strip
                            if '!' in line:
                                idx = line.find('!')
                                cline = line[:idx].strip()
                            elif '*' in line:
                                idx = line.find('*')
                                cline = line[:idx].strip()
                            else:
                                cline.strip()

                            # Cols
                            cols = cline.split()

                            for i,v in enumerate(self.__block2):
                                if v in ivar:
                                    out[v][iz] = float(cols[i])

            # tau?
            if self.__ztype == 1 and 'ltau' in out:
                if out['ltau'][0] <= 0.:
                    out['ltau'][0] = out['ltau'][1]*1e-3
                out['ltau'] = np.log10(out['ltau'])

            # Cut in z?
            if minz is not None:
                ii = np.argmin(out[zvar] - minz)
                for v in out:
                    out[v] = out[v][:ii+1]
            if maxz is not None:
                ii = np.argmin(out[zvar] - maxz)
                for v in out:
                    out[v] = out[v][ii:]

            # If z was not requested
            if needz == 10:
                del out[zvar]

        except:
            raise

        # Return column
        return out

######################################################################
######################################################################
######################################################################

class _atmo_b_1D():
    ''' Class to manage the 1D model atmosphere full format
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_vars': \
          [None,'Get list of variables in the model atmosphere'], \
         'get_vars_units': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding units'], \
         'get_vars_alias': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding alias'], \
          'get_column': \
          [{'minz': \
             'Lower height or optical depth boundary (lower ' + \
             'limit) for output', \
            'maxz': \
             'Upper height or optical depth boundary (upper ' + \
             'limit) for output', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_alias()}'}, \
           'Extract the model atmosphere']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1D model atmosphere
        '''
        try:

            # Number of heights
            self.__nz = 0

            # Open file
            f = open(self.__filename,'r')

            # For each line
            for iline,line in enumerate(f):

                # Strip comments
                if '!' in line:
                    il = line.find('!')
                    if ival < 1: cc = line[il+1:]
                    cline = line[:il]
                elif '*' in line:
                    il = line.find('*')
                    if ival < 1: cc = line[il+1:]
                    cline = line[:il]
                else:
                    cline = line
                    cc = ''

                # Strip spaces
                cline = cline.strip()

                # If empty, skip
                if len(cline) < 1: continue

                # There is content
                self.__nz += 1

            # Initialize variable names
            self.__vars = ['h',r'$\log{\tau_{\rm c}}$', \
                           r'$\chi_{\rm c}$','T',r'P$_{\rm g}$', \
                           r'$\rho$',r'|B|',r'B$_{\rm t}$', \
                           r'B$_{\rm p}$', \
                           r'v$_{\rm x}$',r'v$_{\rm y}$', \
                           r'v$_{\rm z}$', \
                           r'v$_{\rm mi}$',r'P$_{\rm e}$', \
                           r'N$_{\rm e}$',r'N$_{\rm H}$', \
                           r'N$_{\rm H_{\rm a}}$', \
                           r'N$_{\rm H^-}$',r'N$_{\rm H_0}$',
                           r'N$_{\rm H_1}$',r'N$_{\rm H_2}$', \
                           r'N$_{\rm H_3}$',r'N$_{\rm H_4}$', \
                           r'N$_{\rm p^+}$']
            self.__alias = ['h','ltau','chic','T','Pg','rho', \
                            'B','Bt','Bp','vx','vy','vz', \
                            'vmi','Pe','ne','nHT','nHa', \
                            'nH-','nH0','nH1','nH2','nH3','nH4','np']
            self.__vars_units = ['km','','cm^-1]','K','dyn/cm^2', \
                                 'g cm^-3','G','rad','rad','km s^-1', \
                                 'km s^-1','km s^-1','km s^-1',
                                 'dyn cm^-2','cm^-3','cm^-3','cm^-3', \
                                 'cm^-3','cm^-3','cm^-3','cm^-3', \
                                 'cm^-3','cm^-3','cm^-3']
        except:
            raise

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_name(self):
        ''' Get the name in the read file
        '''
        return self.__name

    def _get_comment(self):
        ''' Get the comment in the read file
        '''
        return self.__name

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_vars(self):
        ''' Get variables in model atmosphere
        '''
        return self.__vars

    def _get_wavelength(self):
        ''' Get wavelength for the continuum corresponding to
            the optical depth
        '''
        return self.__wavelength

    def _get_vars_alias(self):
        ''' Get variables and their alias in model atmosphere
        '''
        out = []
        for var,alias in zip(self.__vars,self.__alias):
            out.append(var+' -> '+alias)
        return out

    def _get_vars_units(self):
        ''' Get variables in model atmosphere with units
        '''
        out = []
        for var,uni in zip(self.__vars,self.__vars_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_column(self,minh=None,maxh=None, \
                         mint=None,maxt=None,var=None):
        ''' Get model atmosphere
        '''

        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Initialize variables
        out = {}
        for var in ivar:
            out[var] = np.zeros((self.__nz))

        # Try getting data
        try:

            # Height index
            iz = -1

            # Open file
            with open(self.__filename,'r') as f:

                # Run over lines
                for iline, line in enumerate(f):

                    # Strip
                    if '!' in line:
                        idx = line.find('!')
                        cline = line[:idx].strip()
                    elif '*' in line:
                        idx = line.find('*')
                        cline = line[:idx].strip()
                    else:
                        cline = line.strip()

                    # If empty, skip
                    if len(cline) < 1: continue

                    # Cols
                    cols = cline.split()

                    # Advance
                    iz += 1

                    for i,v in enumerate(self.__alias):
                        if v in ivar:
                            out[v] = col[i]
                        elif v == 'h' or v == 'ltau':
                            if minh is not None or \
                               maxh is not None:
                                out['h'] = col[0]
                            if mint is not None or \
                               maxt is not None:
                                out['tau'] = col[1]

            # Log tau
            if 'ltau' in out:
                if out['ltau'][0] <= 0.:
                    out['ltau'][0] = out['ltau'][1]*1e-3
                out['ltau'] = np.log10(out['ltau'])

            # Cut in z?
            if minh is not None:
                ii = np.argmin(out['h'] - minh)
                for v in out:
                    out[v] = out[v][:ii+1]
            if maxh is not None:
                ii = np.argmin(out['h'] - maxh)
                for v in out:
                    out[v] = out[v][ii:]
            if mint is not None:
                ii = np.argmin(out['h'] - mint)
                for v in out:
                    out[v] = out[v][ii:]
            if maxt is not None:
                ii = np.argmin(out['h'] - maxt)
                for v in out:
                    out[v] = out[v][:ii+1]

            # Remove not requested
            if 'h' in out and 'h' not in ivar:
                del out['h']
            if 'ltau' in out and 'ltau' not in ivar:
                del out['ltau']

        except:
            raise

        # Return column
        return out

######################################################################
######################################################################
######################################################################

class _atmo_15D():
    ''' Class to manage the 3D model atmosphere format
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_precision': \
          [None,'Get type of variable in which the variables ' + \
                'are stored'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height ' + \
                'dimensions'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, and height ' + \
                'dimensions'], \
         'get_vars': \
          [None,'Get list of variables in the model atmosphere'], \
         'get_vars_units': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding units'], \
         'get_vars_alias': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding alias'], \
          'get_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height', \
            'mint': \
             'Lower boundary for output optical depth', \
            'maxt': \
             'Upper boundary for output optical depth', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_alias()}'}, \
           'Extract the model atmosphere for a particular column'], \
          'get_plane': \
          [{'iz': \
             'Coordinate in the height dimension of the ' + \
             'atmospheric parameters to extract', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_alias()}'}, \
           'Extract the model atmosphere for a particular ' + \
           'height index for the whole field of view'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data'] \
          }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 3D model atmosphere
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Read precision
            size = int(struct.unpack('i',f.read(4))[0])
            if size == 4:
                self.__vtype = np.float32
                self.__fmt = 'f'
                self.__byt = 4
            elif size == 8:
                self.__vtype = np.float64
                self.__fmt = 'd'
                self.__byt = 8
            else:
                _error(f"Size {size} not valid",1)
                return False

            # Get dimensions
            self.__nx = int(struct.unpack('i',f.read(4))[0])
            self.__ny = int(struct.unpack('i',f.read(4))[0])
            self.__nz = int(struct.unpack('i',f.read(4))[0])

            # Close file
            f.close()

            # Size of head
            self.__head = 4 + 4 + 4*3

            # Variables
            self.__nvar = 24
            self.__vars = ['h',r'$\log{\tau_{\rm c}}$', \
                           r'$\chi_{\rm c}$', \
                           'T',r'P$_{\rm g}$',r'$\rho$', \
                           r'B$_{\rm x}$',r'B$_{\rm y}$', \
                           r'B$_{\rm z}$', \
                           r'v$_{\rm x}$',r'v$_{\rm y}$', \
                           r'v$_{\rm z}$', \
                           r'v$_{\rm mi}$',r'P$_{\rm e}$', \
                           r'N$_{\rm e}$',r'N$_{\rm H}$', \
                           r'N$_{\rm H_{\rm a}}$', \
                           r'N$_{\rm H^-}$',r'N$_{\rm H_0}$',
                           r'N$_{\rm H_1}$',r'N$_{\rm H_2}$', \
                           r'N$_{\rm H_3}$',r'N$_{\rm H_4}$', \
                           r'N$_{\rm p^+}$']
            self.__alias = ['h','ltau','chic','T','Pg','rho', \
                            'Bx','By','Bz','vx','vy','vz', \
                            'vmi','Pe','ne','nHT','nHa', \
                            'nH-','nH0','nH1','nH2','nH3','nH4','np']
            self.__vars_units = ['[km]','',r'[cm$^{-1}$]', \
                                 '[K]',r'[dyn cm$^{-2}$]', \
                                r'[g cm$^{-3}$]','[G]','[G]','[G]', \
                                r'[km s$^{-1}$]',r'[km s$^{-1}$]', \
                                r'[km s$^{-1}$]',r'[km s$^{-1}$]', \
                                r'[dyn cm$^{-2}$]',r'[cm$^{-3}$]', \
                                r'[cm$^{-3}$]',r'[cm$^{-3}$]', \
                                r'[cm$^{-3}$]',r'[cm$^{-3}$]', \
                                r'[cm$^{-3}$]',r'[cm$^{-3}$]', \
                                r'[cm$^{-3}$]',r'[cm$^{-3}$]', \
                                r'[cm$^{-3}$]']
            # Size of column
            self.__column = self.__byt*self.__nz*self.__nvar

            # Create memmaps
            self.__data = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__head, \
                                    dtype=self.__vtype, \
                                    shape=(self.__nx, \
                                           self.__ny, \
                                           self.__nvar, \
                                           self.__nz))

        except struct.error:
            return False
            raise
        except:
            raise

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_precision(self):
        ''' Get type of variable in which the variables are stored
        '''
        return self.__vtype

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nxyz(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self.__nx, self.__ny, self.__nz

    def _get_dims(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self._get_nxyz()

    def _get_vars(self):
        ''' Get variables in model atmosphere
        '''
        return self.__vars

    def _get_vars_alias(self):
        ''' Get variables and their alias in model atmosphere
        '''
        out = []
        for var,alias in zip(self.__vars,self.__alias):
            out.append(var+' -> '+alias)
        return out

    def _get_vars_units(self):
        ''' Get variables in model atmosphere with units
        '''
        out = []
        for var,uni in zip(self.__vars,self.__vars_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_column(self,ix,iy,minh=None,maxh=None, \
                         mint=None,maxt=None,var=None):
        ''' Get model atmosphere for a given pixel
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Try getting data
        try:

            # Read data
            col = self.__data[ix,iy,:,:]

            # Adjust height
            if minh is not None:
                i = np.argmin(np.absolute(col[0,:] - minh))
                col = col[:,:i+1]
            if maxh is not None:
                i = np.argmin(np.absolute(col[0,:] - maxh))
                col = col[:,i:]
            if mint is not None:
                i = np.argmin(np.absolute(col[1,:] - mint))
                col = col[:,i:]
            if maxt is not None:
                i = np.argmin(np.absolute(col[1,:] - maxt))
                col = col[:,:i+1]

        except:
            raise

        # Return column
        out = {}
        for i,v in enumerate(self.__alias):
            if v in ivar:
                out[v] = col[i,:].copy()
                if v == 'ltau':
                    if out[v][0] <= 0:
                        out[v][0] = out[v][1]*1e-3
                    out[v] = np.log10(out[v])
        return out

    def _get_plane(self,iz,var=None):
        ''' Get model atmosphere for a given height index
        '''

        # Valid?
        if not isinstance(iz, int) and not isinstance(iz, npint):
           _error('iz must be an integer',1)
           return None
        if iz < 0 or iz >= self.__nz:
           _error('The requested height index is out of bounds',1)
           return None


        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Try getting data
        try:

            # Get plane
            col = self.__data[:,:,:,iz]

        except:
            raise

        # Return plane
        out = {}
        for i,v in enumerate(self.__alias):
            if v in ivar:
                out[v] = col[:,:,i].copy()
                if v == 'ltau':
                    out[v] = np.where(out[v] <= 0., \
                                      1e-16,out[v])
                    out[v] = np.log10(out[v])

        return out

    def _get_cube(self):
        ''' Get whole model
        '''
        return self.__data

######################################################################
######################################################################
######################################################################

class _inversion_in():
    ''' Class to manage the input data file for the TIC
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_info_verb': \
          [None,'Get verbosity information about the file'], \
         'get_info': \
          [None,'Get information integers about the file'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, and ' + \
                'wavelength dimensions'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
         'get_los': \
          [None,'Get the line of sight'], \
         'get_los_column': \
          [None,'Get the line of sight for a given pixel'], \
         'get_los_plane': \
          [None,'Get the line of sight for the whole plane'], \
         'get_stokesi_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the intensity at a particular column'], \
         'get_stokesq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to ' + \
             'extract', \
            'iy': \
             'Coordinate in the y dimension of the column to ' + \
             'extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes Q at a particular column'], \
         'get_stokesu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes U at a particular column'], \
         'get_stokesv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes V at a particular column'], \
         'get_linear_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract linear polarization at a particular column'], \
         'get_stokes_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to ' + \
             'extract', \
            'iy': \
             'Coordinate in the y dimension of the column to ' + \
             'extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes parameters at a particular column'], \
         'get_stokesi_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract Stokes I at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_linear_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a ' + \
           'particular wavelength index for the whole ' + \
           'field of view'], \
         'get_stokes_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular ' + \
           'wavelength index for ' + \
           'the whole field of view'], \
         'get_sigma': \
          [None,'Get the sigma for Stokes parameters'], \
         'get_sigmai_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the intensity sigma at a particular column'], \
         'get_sigmaq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to ' + \
             'extract', \
            'iy': \
             'Coordinate in the y dimension of the column to ' + \
             'extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes Q sigma at a particular column'], \
         'get_sigmau_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes U sigma at a particular column'], \
         'get_sigmav_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes V sigma at a particular column'], \
         'get_sigma_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes sigma at a particular column'], \
         'get_sigmai_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the intensity sigma at a particular ' + \
           'wavelength index for ' + \
           'the whole field of view'], \
         'get_sigmaq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes Q sigma at a particular ' + \
           'wavelength index for ' + \
           'the whole field of view'], \
         'get_sigmau_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes U sigma at a particular ' + \
           'wavelength index for ' + \
           'the whole field of view'], \
         'get_sigmav_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes V sigma at a particular ' + \
           'wavelength index for ' + \
           'the whole field of view'], \
         'get_sigma_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes sigma at a particular ' + \
           'wavelength index for ' + \
           'the whole field of view'], \
         'get_diff': \
          [None,'Get the diffuse light Stokes parameters'], \
         'get_diffi_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the intensity diffuse light at a ' + \
           'particular column'], \
         'get_diffq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes Q diffuse light at a ' + \
           'particular column'], \
         'get_diffu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes U diffuse light at a particular ' + \
           'column'], \
         'get_diffv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes V diffuse light at a particular ' + \
           'column'], \
         'get_diff_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes diffuse light at a ' + \
           'particular column'], \
         'get_diffi_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the intensity diffuse light at a ' + \
           'particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diffq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes Q diffuse light at a ' + \
           'particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diffu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes U diffuse light at a ' + \
           'particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diffv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes V diffuse light at a ' + \
           'particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diff_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes diffuse light at a ' + \
           'particular wavelength index for ' + \
           'the whole field of view'], \
          'get_cube': \
          [None,f'Get a memmap to the pixel wise data.'] \
          }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert inversion output file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Dimensions
            self.__nx = struct.unpack('i',f.read(4))[0]
            self.__ny = struct.unpack('i',f.read(4))[0]
            self.__nl = struct.unpack('i',f.read(4))[0]

            # Read info
            self.__info = np.array(struct.unpack('iiii',f.read(16)))

            # Close
            f.close()

        except:
            raise

        # Size header
        self.__to_lambda = 4*8
        self.__head = self.__to_lambda + self.__nl*8
        self.__to_data = self.__head

        # Sanity
        if self.__info[0] < 0 or self.__info[0] > 1:
            _error('Stokes label must be 0 or 1',1)
            return False
        if self.__info[1] < 0 or self.__info[1] > 1:
            _error('LOS label must be 0 or 1',1)
            return False
        if self.__info[2] < 0 or self.__info[2] > 4:
            _error('Sigma label must be 0, 1, 2, 3, or 4',1)
            return False
        if self.__info[3] < 0 or self.__info[3] > 4:
            _error('Diffuse light label must be 0, 1, 2, 3, or 4',1)
            return False

        # If constant LOS
        if self.__info[1] == 0:
            self.__to_los = self.__to_data
            self.__to_data += 16
        # Variable LOS
        elif self.__info[1] == 1:
            self.__to_los = 0

        # No sigma
        if self.__info[2] == 0:
            self.__to_sigma = self.__to_data
        # If constant sigma intensity
        elif self.__info[2] == 1 and self.__info[0] == 0:
            self.__to_sigma = self.__to_data
            self.__to_data += 8
        # If constant sigma polarization
        elif self.__info[2] == 1 and self.__info[0] == 1:
            self.__to_sigma = self.__to_data
            self.__to_data += 8*4
        # If variable sigma intensity
        elif self.__info[2] == 2 and self.__info[0] == 0:
            self.__to_sigma = self.__to_data
            self.__to_data += 8*self.__nl
        # If variable sigma polarization
        elif self.__info[2] == 2 and self.__info[0] == 1:
            self.__to_sigma = self.__to_data
            self.__to_data += 8*self.__nl*4
        # If pixel constant sigma intensity
        elif self.__info[2] == 3 and self.__info[0] == 0:
            self.__to_sigma = 0
        # If pixel constant sigma polarization
        elif self.__info[2] == 3 and self.__info[0] == 1:
            self.__to_sigma = 0
        # If pixel variable sigma intensity
        elif self.__info[2] == 4 and self.__info[0] == 0:
            self.__to_sigma = 0
        # If pixel variable sigma polarization
        elif self.__info[2] == 4 and self.__info[0] == 1:
            self.__to_sigma = 0

        # No diffuse light
        if self.__info[3] == 0:
            self.__to_diff = self.__to_data
        # If constant diffuse light intensity
        elif self.__info[3] == 1:
            self.__to_diff = self.__to_data
            self.__to_data += 8*self.__nl
        # If constant diffuse light polarization
        elif self.__info[3] == 2:
            self.__to_diff = self.__to_data
            self.__to_data += 8*4
        # If pixel diffuse light intensity
        elif self.__info[3] == 3:
            self.__to_diff = 0
        # If pixel diffuse light polarization
        elif self.__info[3] == 4:
            self.__to_diff = 0


        #
        # Build type for the data block
        #

        # Initialize
        fields = []

        # LOS
        if self.__info[1] == 1:
            fields.append(('los',np.float64,(2)))

        # Stokes
        # Only intensity
        if self.__info[0] == 0:
            ns = 1
        # Polarization
        else:
            ns = 4
        fields.append(('stokes',np.float64,(ns,self.__nl)))

        # Sigma
        # Constant
        if self.__info[2] == 3:
            fields.append(('sigma',np.float64,(ns,1)))
        # Profile
        elif self.__info[2] == 4:
            fields.append(('sigma',np.float64,(ns,self.__nl)))

        # Diffuse light
        # Intensity
        if self.__info[3] == 3:
            fields.append(('diff',np.float64,(1,self.__nl)))
        # Profile
        elif self.__info[3] == 4:
            fields.append(('diff',np.float64,(4,self.__nl)))

        # Data type
        self.__data_dtype = np.dtype(fields)


        #
        # Create memmaps
        #

        # Wavelength
        self.__lam = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        # Constant LOS
        if self.__info[1] == 0:
            self.__los_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_los, \
                                      dtype=np.float64, \
                                      shape=(2))

        # If constant sigma intensity
        if self.__info[2] == 1 and self.__info[0] == 0:
            self.__sig_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_sigma, \
                                      dtype=np.float64, \
                                      shape=(1))
        # If constant sigma polarization
        elif self.__info[2] == 1 and self.__info[0] == 1:
            self.__sig_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_sigma, \
                                      dtype=np.float64, \
                                      shape=(4))
        # If constant but profile sigma intensity
        elif self.__info[2] == 2 and self.__info[0] == 0:
            self.__sig_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_sigma, \
                                      dtype=np.float64, \
                                      shape=(self.__nl))
        # If constant but profile sigma polarization
        elif self.__info[2] == 1 and self.__info[0] == 1:
            self.__sig_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_sigma, \
                                      dtype=np.float64, \
                                      shape=(4,self.__nl))

        # If constant diffuse light intensity
        if self.__info[3] == 1:
            self.__dif_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_diff, \
                                      dtype=np.float64, \
                                      shape=(self.__nl))
        # If constant diffuse light polarization
        elif self.__info[3] == 2:
            self.__dif_ct = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__to_diff, \
                                      dtype=np.float64, \
                                      shape=(4,self.__nl))

        # Data
        self.__data = np.memmap(self.__filename, \
                                mode='r', \
                                offset=self.__to_data, \
                                dtype=self.__data_dtype, \
                                shape=(self.__nx,self.__ny))

        # Return valid
        return True


    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_info_verb(self):
        ''' Get verbosity information about the file
        '''
        if self.__info[0] == 0:
            _verbose('Only intensity file',0)
        else:
            _verbose('Full Stokes file',0)
        if self.__info[1] == 0:
            _verbose('Single LOS',0)
        else:
            _verbose('Pixelwise LOS',0)
        if self.__info[2] == 0:
            _verbose('No sigma',0)
        elif self.__info[2] == 1:
            _verbose('Constant sigma',0)
        elif self.__info[2] == 2:
            _verbose('Wavelength dependent constant sigma',0)
        elif self.__info[2] == 3:
            _verbose('Constant pixelwise sigma',0)
        elif self.__info[2] == 4:
            _verbose('Wavelength dependent pixelwise sigma',0)
        if self.__info[3] == 0:
            _verbose('No diffuse light profile',0)
        elif self.__info[3] == 1:
            _verbose('Constant only intensity diffuse ' + \
                     'light profile',0)
        elif self.__info[3] == 2:
            _verbose('Constant full Stokes diffuse light profile',0)
        elif self.__info[3] == 3:
            _verbose('Pixelwise only intensity diffuse ' + \
                     'light profile',0)
        elif self.__info[3] == 4:
            _verbose('Pixelwise full Stokes diffuse light profile',0)

    def _get_info(self,i=None):
        ''' Get index in information
        '''
        if i is None:
            return self.__info
        else:
            return self.__info[i]

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in x, y, and wavelength axes
        '''
        return self.__nx, self.__ny, self.__nl

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = self.__lam.copy()
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def __get_los_ct(self):
        ''' Get constant los
        '''
        try:
            los = self.__los_ct.copy()
        except:
            raise
        return los

    def __get_los_column(self,ix,iy):
        ''' Get los at a column
        '''
        try:
            los = self.__data['los'][ix,iy].copy()
        except:
            raise
        return los

    def __get_los_plane(self):
        ''' Get los for the plane
        '''
        try:
            los = self.__data['los'].copy()
        except:
            raise
        return los

    def _get_los(self):
        ''' Get LOS if full constant
        '''
        # Constant
        if self.__info[1] == 0:
            return self.__get_los_ct()
        else:
            return self.__get_los_plane()

    def _get_los_column(self,ix=None,iy=None):
        ''' Get LOS at a given column
        '''
        # Constant
        if self.__info[1] == 0:
            return self.__get_los_ct()
        # Variable
        else:
            # Valid?
            if not isinstance(ix, int) and not isinstance(ix, npint):
               _error('ix must be an integer',1)
               return None
            if not isinstance(iy, int) and not isinstance(iy, npint):
               _error('iy must be an integer',1)
               return None
            if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
               _error('The requested column is out of bounds',1)
               return None
            return self.__get_los_column(ix,iy)

    def _get_los_plane(self):
        ''' Get LOS for the whole FoV
        '''
        # Constant
        if self.__info[1] == 0:
            los = np.empty((self.__nx,self.__ny,2))
            llos = self.__get_los_ct()
            los[:,:,0] = llos[0]
            los[:,:,1] = llos[1]
            return los
        # Variable
        else:
            return self.__get_los_plane()

    def __get_gen_column(self,ix,iy,field,pol, \
                         minl=None,maxl=None,nl=None, \
                         fractional=False,indx=[0]):
        ''' Generic read of a column parameter
        '''

        # Output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Try getting data
        try:

            # Intensity
            if 0 in indx or fractional:

                # Get intensity
                arrI = self.__data[field][ix,iy,0,:].copy()

                # Out?
                if 0 in indx:
                    out[0] = arrI

            # Q, U, and V
            for j in range(1,4):

                # There is polarization
                if pol:

                    # To output
                    if j in indx:

                        # Read Stokes
                        out[j] = self.__data[field][ix,iy,j,:].copy()

                # No pol
                else:

                    # Zeros
                    out[j] = np.zeros((nl))

                # Manage units
                if fractional and j in indx:
                    out[j] /= arrI

            # Adjust wavelength
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                for j in indx:
                    out[j] = out[j][i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                for j in indx:
                    out[j] = out[j][:i+1]

        # Others
        except:
            raise

        # Return
        return out

    def _get_stokesi_column(self,ix,iy,minl=None,maxl=None):
        ''' Get intensity profile at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'stokes', \
                                     self.__info[0]==1, \
                                     minl,maxl,self.__nl, \
                                     False,[0])[0]

    def _get_stokesq_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes Q profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
           return np.zeros((self.__nl))

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'stokes', \
                                     self.__info[0]==1, \
                                     minl,maxl,self.__nl, \
                                     fractional,[1])[1]

    def _get_stokesu_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes U profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
           return np.zeros((self.__nl))

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'stokes', \
                                     self.__info[0]==1, \
                                     minl,maxl,self.__nl, \
                                     fractional,[2])[2]

    def _get_stokesv_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes V profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
           return np.zeros((self.__nl))

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'stokes', \
                                     self.__info[0]==1, \
                                     minl,maxl,self.__nl, \
                                     fractional,[3])[3]

    def _get_linear_column(self,ix,iy,minl=None,maxl=None, \
                           fractional=False):
        ''' Get Stokes linear polarization profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
           return np.zeros((self.__nl))

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        qu = self.__get_gen_column(ix,iy,'stokes', \
                                   self.__info[0]==1, \
                                   minl,maxl,self.__nl, \
                                   fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_column(self,ix,iy,minl=None,maxl=None, \
                           fractional=False):
        ''' Get Stokes parameter at a given column
        '''


        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # Mode?
        if self.__info[0] == 0:
            lout = self._get_stokesi_column(ix,iy, \
                                            minl=minl, \
                                            maxl=maxl)
            nl = lout.size
            out = np.zeros((4,nl))
            out[0] = lout
            return out
        else:
            iquv = self.__get_gen_column(ix,iy,'stokes', \
                                         self.__info[0]==1, \
                                         minl,maxl,self.__nl, \
                                         fractional,[0,1,2,3])
            return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))


    def __get_gen_plane(self,il,field,pol, \
                        fractional=False,indx=[0]):
        ''' Generic read of a plane parameter
        '''

        # Output
        out = [None,None,None,None]

        # Try getting data
        try:

            # Intensity
            if 0 in indx or fractional:

                # Get intensity
                arrI = self.__data[field][:,:,0,il].copy()

                # Out?
                if 0 in indx:
                    out[0] = arrI

            # Q, U, and V
            for j in range(1,4):

                # There is polarization
                if pol:

                    # To output
                    if j in indx:

                        # Get Stokes
                        out[j] = self.__data[field][:,:,j,il].copy()

                # No polarization
                else:

                    # Zero
                    if j in indx:
                        out[j] = np.zeros((self.__nx,self.__ny))

                # Manage units
                if fractional and j in indx:
                    out[j] /= arrI
        except:
            raise

        # Return
        return out

    def _get_stokesi_plane(self,il):
        ''' Get intensity profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'stokes', \
                                    self.__info[0]==1, \
                                    False,[0])[0]

    def _get_stokesq_plane(self,il,fractional=False):
        ''' Get Stokes Q profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'stokes', \
                                    self.__info[0]==1, \
                                    fractional,[1])[1]

    def _get_stokesu_plane(self,il,fractional=False):
        ''' Get Stokes U profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'stokes', \
                                    self.__info[0]==1, \
                                    fractional,[2])[2]

    def _get_stokesv_plane(self,il,fractional=False):
        ''' Get Stokes V profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'stokes', \
                                    self.__info[0]==1, \
                                    fractional,[3])[3]

    def _get_linear_plane(self,il,fractional=False):
        ''' Get Stokes V profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        qu = self.__get_gen_plane(il,'stokes', \
                                  self.__info[0]==1, \
                                  fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_plane(self,il,fractional=False):
        ''' Get Stokes profiles at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        iquv = self.__get_gen_plane(il,'stokes', \
                                    self.__info[0]==1, \
                                    fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def __no_sigma(self):
        ''' Message that there is no sigma
        '''
        _error('There is no sigma in this file',1)
        return None

    def __get_sigma_ct(self):
        ''' Get constant sigma
        '''
        try:
            # Full constant
            if self.__info[2] == 1:
                if self.__info[0] == 0:
                    sig = np.zeros((4))
                    sig[0] = self.__sig_ct.copy()[0]
                else:
                    sig = self.__sig_ct.copy()
            # Wavelength dependent constant
            elif self.__info[2] == 2:
                if self.__info[0] == 0:
                    sig = np.zeros((4,self.__nl))
                    sig[0,:] = self.__sig_ct.copy()
                else:
                    sig = self.__sig_ct.copy()
        except:
            raise
        return sig

    def _get_sigma(self):
        ''' Get sigma if full constant
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant
        if self.__info[2] == 1 or self.__info[2] == 2:
            return self.__get_sigma_ct()
        # Constant pixelwise
        elif self.__info[2] == 3:
            iquv = self.__get_gen_plane(0,'sigma', \
                                        self.__info[0]==1, \
                                        False,[0,1,2,3])
            return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))
        # Other
        else:
            _error('Sigma is pixelwise and non-constant, ' + \
                   'get a column with get_sigma_column() ' + \
                   'or a plane with get_sigma_plane()',1)
            return None

    def _get_sigmai_column(self,ix,iy,minl=None,maxl=None):
        ''' Get sigma for the intensity at a given column
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1 or self.__info[2] == 2:
            return self.__get_sigma_ct()[0]
        if self.__info[2] == 3:
            nl = 1
            iminl = None
            imaxl = None
        elif self.__info[2] == 4:
            nl = self.__nl
            iminl = minl
            imaxl = maxl

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'sigma', \
                                     self.__info[0]==1, \
                                     iminl,imaxl,nl,False,[0])[0]

    def _get_sigmaq_column(self,ix,iy,minl=None,maxl=None):
        ''' Get sigma for Stokes Q at a given column
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1 or self.__info[2] == 2:
            return self.__get_sigma_ct()[1]
        if self.__info[2] == 3:
            nl = 1
            iminl = None
            imaxl = None
        elif self.__info[2] == 4:
            nl = self.__nl
            iminl = minl
            imaxl = maxl

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'sigma', \
                                     self.__info[0]==1, \
                                     iminl,imaxl,nl,False,[1])[1]

    def _get_sigmau_column(self,ix,iy,minl=None,maxl=None):
        ''' Get sigma for Stokes U at a given column
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1 or self.__info[2] == 2:
            return self.__get_sigma_ct()[2]
        if self.__info[2] == 3:
            nl = 1
            iminl = None
            imaxl = None
        elif self.__info[2] == 4:
            nl = self.__nl
            iminl = minl
            imaxl = maxl

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'sigma', \
                                     self.__info[0]==1, \
                                     iminl,imaxl,nl,False,[2])[2]

    def _get_sigmav_column(self,ix,iy,minl=None,maxl=None):
        ''' Get sigma for Stokes V at a given column
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1 or self.__info[2] == 2:
            return self.__get_sigma_ct()[3]
        if self.__info[2] == 3:
            nl = 1
            iminl = None
            imaxl = None
        elif self.__info[2] == 4:
            nl = self.__nl
            iminl = minl
            imaxl = maxl

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'sigma', \
                                     self.__info[0]==1, \
                                     iminl,imaxl,nl,False,[3])[3]

    def _get_sigma_column(self,ix,iy,minl=None,maxl=None):
        ''' Get sigma for the intensity at a given column
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1 or self.__info[2] == 2:
            return self.__get_sigma_ct()
        if self.__info[2] == 3:
            nl = 1
            iminl = None
            imaxl = None
        elif self.__info[2] == 4:
            nl = self.__nl
            iminl = minl
            imaxl = maxl

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        iquv = self.__get_gen_column(ix,iy,'sigma', \
                                     self.__info[0]==1, \
                                     iminl,imaxl,nl,False,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def _get_sigmai_plane(self,il):
        ''' Get sigma for the intensity at a given wavelength
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[0]
            return sig
        elif self.__info[2] == 2:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[0][il]
            return sig
        if self.__info[2] == 3:
            nl = 1
        elif self.__info[2] == 4:
            nl = self.__nl

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'sigma', \
                                    self.__info[0]==1, \
                                    False,[0])[0]

    def _get_sigmaq_plane(self,il):
        ''' Get sigma for Stokes Q at a given wavelength
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[1]
            return sig
        elif self.__info[2] == 2:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[1][il]
            return sig
        if self.__info[2] == 3:
            nl = 1
        elif self.__info[2] == 4:
            nl = self.__nl

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'sigma', \
                                    self.__info[0]==1, \
                                    False,[1])[1]

    def _get_sigmau_plane(self,il):
        ''' Get sigma for Stokes U at a given wavelength
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[2]
            return sig
        elif self.__info[2] == 2:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[2][il]
            return sig
        if self.__info[2] == 3:
            nl = 1
        elif self.__info[2] == 4:
            nl = self.__nl

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'sigma', \
                                    self.__info[0]==1, \
                                    False,[2])[2]

    def _get_sigmav_plane(self,il):
        ''' Get sigma for Stokes V at a given wavelength
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[3]
            return sig
        elif self.__info[2] == 2:
            sig = np.empty((self.__nx,self.__ny))
            sig[:,:] = self.__get_sigma_ct()[3][il]
            return sig
        if self.__info[2] == 3:
            nl = 1
        elif self.__info[2] == 4:
            nl = self.__nl

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'sigma', \
                                    self.__info[0]==1, \
                                    False,[3])[3]

    def _get_sigma_plane(self,il):
        ''' Get sigma at a given wavelength
        '''
        # No sigma
        if self.__info[2] == 0: return self.__no_sigma()
        # Constant with pixel
        if self.__info[2] == 1:
            sig = np.empty((4,self.__nx,self.__ny))
            lsig = self.__get_sigma_ct()
            for i in range(4):
                sig[i,:,:] = lsig[i]
            return sig
        elif self.__info[2] == 2:
            sig = np.empty((4,self.__nx,self.__ny))
            lsig = self.__get_sigma_ct()
            for i in range(4):
                sig[i,:,:] = lsig[i][il]
            return sig
        if self.__info[2] == 3:
            nl = 1
        elif self.__info[2] == 4:
            nl = self.__nl

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        iquv = self.__get_gen_plane(il,'sigma', \
                                    self.__info[0]==1, \
                                    False,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def __no_diff(self):
        ''' Message that there is no diffuse light
        '''
        _error('There is no diffuse light in this file',1)
        return None

    def __get_diff_ct(self):
        ''' Get constant diffuse light
        '''
        try:
            # Intensity
            if self.__info[3] == 1:
                diff = np.zeros((4,self.__nl))
                diff[0,:] = self.__dif_ct.copy()
            # Polarized
            elif self.__info[3] == 2:
                diff = self.__dif_ct.copy()
        except:
            raise
        return diff

    def _get_diff(self):
        ''' Get diff if full constant
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant
        if self.__info[3] == 1 or self.__info[3] == 2:
            return self.__get_diff_ct()
        # Other
        else:
            _error('Diffuse light is pixelwise, get a ' + \
                   'column with get_diff_column() or ' + \
                   'a plane with get_diff_plane()',1)
            return None

    def _get_diffi_column(self,ix,iy,minl=None,maxl=None):
        ''' Get diffuse light for the intensity at a given column
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            return self.__get_diff_ct()[0]

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'diff', \
                                     self.__info[3]==4, \
                                     minl,maxl,self.__nl, \
                                     False,[0])[0]

    def _get_diffq_column(self,ix,iy,minl=None,maxl=None):
        ''' Get diff for Stokes Q at a given column
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            return self.__get_diff_ct()[1]

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'diff', \
                                     self.__info[3]==4, \
                                     minl,maxl,self.__nl, \
                                     False,[1])[1]

    def _get_diffu_column(self,ix,iy,minl=None,maxl=None):
        ''' Get diff for Stokes U at a given column
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            return self.__get_diff_ct()[2]

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'diff', \
                                     self.__info[3]==4, \
                                     minl,maxl,self.__nl, \
                                     False,[2])[2]

    def _get_diffv_column(self,ix,iy,minl=None,maxl=None):
        ''' Get diff for Stokes V at a given column
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            return self.__get_diff_ct()[3]

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,'diff', \
                                     self.__info[3]==4, \
                                     minl,maxl,self.__nl, \
                                     False,[3])[3]

    def _get_diff_column(self,ix,iy,minl=None,maxl=None):
        ''' Get diff for the intensity at a given column
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            return self.__get_diff_ct()

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        iquv = self.__get_gen_column(ix,iy,'diff', \
                                     self.__info[3]==4, \
                                     minl,maxl,self.__nl, \
                                     False,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def _get_diffi_plane(self,il):
        ''' Get diffuse light for the intensity at a given wavelength
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            diff = np.empty((self.__nx,self.__ny))
            diff[:,:] = self.__get_diff_ct()[0][il]

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'diff', \
                                    self.__info[3]==4, \
                                    False,[0])[0]

    def _get_diffq_plane(self,il):
        ''' Get diffuse light for Stokes Q at a given wavelength
        '''
        # No diff
        if self.__info[2] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            diff = np.empty((self.__nx,self.__ny))
            diff[:,:] = self.__get_diff_ct()[1][il]

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'diff', \
                                    self.__info[3]==4, \
                                    False,[1])[1]

    def _get_diffu_plane(self,il):
        ''' Get diffuse light for Stokes U at a given wavelength
        '''
        # No diff
        if self.__info[2] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            diff = np.empty((self.__nx,self.__ny))
            diff[:,:] = self.__get_diff_ct()[2][il]

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'diff', \
                                    self.__info[3]==4, \
                                    False,[2])[2]

    def _get_diffv_plane(self,il):
        ''' Get diffuse light for Stokes V at a given wavelength
        '''
        # No diff
        if self.__info[2] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            diff = np.empty((self.__nx,self.__ny))
            diff[:,:] = self.__get_diff_ct()[3][il]

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        return self.__get_gen_plane(il,'diff', \
                                    self.__info[3]==4, \
                                    False,[3])[3]

    def _get_diff_plane(self,il):
        ''' Get diffuse light for full Stokes at a given wavelength
        '''
        # No diff
        if self.__info[2] == 0: return self.__no_diff()
        # Constant with pixel
        if self.__info[3] == 1 or self.__info[3] == 2:
            diff = np.empty((4,self.__nx,self.__ny))
            ldiff = self.__get_diff_ct()
            for i in range(4):
                diff[i,:,:] = ldiff[i][il]

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        iquv = self.__get_gen_plane(il,'diff', \
                                    self.__info[3]==4, \
                                    False,[0,1,2,3])

        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def _get_cube(self):
        ''' Get the whole data
        '''
        return self.__data

######################################################################
######################################################################
######################################################################

class _inversion_out():
    ''' Class to manage the Result file from the TIC
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_polarized': \
          [None,'Get if the inversion included polarization'], \
         'get_vtype': \
          [None,'Get the type of reference frame for the ' + \
                'velocity vector'], \
         'get_btype': \
          [None,'Get the type of reference frame for the ' +  \
                'magnetic field vector'], \
         'get_jkqin': \
          [None,'Get if there are inverted ad-hoc JKQ'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and ' + \
                'height dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_nvar_atmo': \
          [None,'Get number of variables in model atmosphere'], \
         'get_nvar': \
          [None,'Get number of variables in node data'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, height, ' + \
                'and wavelength dimensions'], \
         'get_vars': \
          [None,'Get list of variables with available ' + \
                'node results'], \
         'get_vars_units': \
          [None,'Get list of variables with available node ' + \
                'results with their corresponding units'], \
         'get_vars_atmo': \
          [None,'Get list of variables in the model atmosphere'], \
         'get_vars_atmo_units': \
          [None,'Get list of variables in the model atmosphere ' + \
                'with their corresponding units'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars()}'}, \
           'Extract the result of the inversion at a particular ' + \
           'column'], \
          'get_column_atmo': \
          [{'ix': \
             'Coordinate in the x dimension of the column to ' + \
             'extract', \
            'iy': \
             'Coordinate in the y dimension of the column to ' + \
             'extract', \
            'minh': \
             'Lower boundary for output optical depth', \
            'maxh': \
             'Upper boundary for output optical depth', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_atmo()}'}, \
           'Extract the model atmosphere for a particular column'], \
          'get_column_rf': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract the response function at a particular column'], \
         'get_plane_chi': \
          [None,'Get the value of the initial and final merit ' + \
                'function for the whole field of view'], \
          'get_plane_stk': \
          [{'il': \
             'Coordinate in the wavelength dimension of the ' + \
             'Stokes parameters to extract', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract Stokes parameters at a given wavelength ' + \
           'position for the whole field of view'], \
          'get_plane_atmo': \
          [{'iz': \
             'Coordinate in the height dimension of the ' + \
             'atmospheric parameters to extract', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_atmo()}'}, \
           'Extract the model atmosphere for a particular ' + \
           'height index for the whole field of view'], \
          'get_node': \
          [{'var': \
             'Variables for which to extract the node information'}, \
           'Extract the node information (full cube) for a ' + \
           'given variable'], \
          'get_cube': \
          [None,f'Get a memmap to the fit data. The file must ' + \
           'be complete to use this method.'], \
          'get_cube_atmo': \
          [None,f'Get a memmap to the whole inverted model ' + \
           'atmosphere. The file must be complete to use ' + \
           'this method.']
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self,debug=False):
        ''' Reads hanlert inversion output file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Read info inversion
            info = struct.unpack('i',f.read(4))[0]
            if debug: print(f'info {info}')

            # JKQ
            self.__jkqa = info > 7
            if debug: print(f'JKQ info {self.__jkqa}')
            if info > 7: info -= 8

            # v type
            if info > 3:
                self.__vtype = 1
                info -= 4
                self.__vx = 'vpos'
                self.__vy = 'vlphi'
                self.__vz = 'vlos'
                self.__vxu = 'km s^-1'
                self.__vyu = 'rad'
                self.__vzu = 'km s^-1'
            else:
                self.__vtype = 0
                self.__vx = 'vx'
                self.__vy = 'vy'
                self.__vz = 'vz'
                self.__vxu = 'km s^-1'
                self.__vyu = 'km s^-1'
                self.__vzu = 'km s^-1'

            # B type
            if info > 1:
                self.__btype = 1
                info -= 2
                self.__bx = 'Blos'
                self.__by = 'Btra'
                self.__bz = 'Blphi'
                self.__bxu = 'G'
                self.__byu = 'G'
                self.__bzu = 'rad'
            else:
                self.__btype = 0
                self.__bx = 'B'
                self.__by = 'Btheta'
                self.__bz = 'Bphi'
                self.__bxu = 'G'
                self.__byu = 'rad'
                self.__bzu = 'rad'

            # Polarization
            self.__polarization = info > 0
            if debug: print(f'Polarization {self.__polarization}')
            if debug: print(f'Btype {self.__btype}')
            if debug: print(f'vtype {self.__vtype}')

            # Dimensions
            self.__nx = struct.unpack('i',f.read(4))[0]
            self.__ny = struct.unpack('i',f.read(4))[0]
            self.__nz = struct.unpack('i',f.read(4))[0]
            if debug: print(f'nx {self.__nx}')
            if debug: print(f'ny {self.__ny}')
            if debug: print(f'nz {self.__nz}')

            #
            # Prepare jumps atmosphere related

            # Header
            self.__s_head = 5*4

            # Jump to beginning of atmosphere
            self.__jump_to_atmo = self.__s_head

            # Number of variables in atmospheric model
            if self.__jkqa:
              self.__nvar_atmo = 27
              self.__vars_atmo = ['ltau','T','Pg', \
                                  'Bx','By','Bz','vx', \
                                  'vy','vz','vturb','ne', \
                                  'nHT','nHa','nH(0)','nH(1)', \
                                  'nH(2)','nH(3)','nH(4)','np', \
                                  'J10','Re{J11}','Im{J11}', \
                                  'J20','Re{J21}','Im{J21}', \
                                  'Re{J22}','Im{J22}','f']
              self.__vars_atmo_units = ['','K','dyn cm^-2','G', \
                                        'G','G','km s^-1','km s^-1', \
                                        'km s^-1','km s^-1','cm^-3', \
                                        'cm^-3','cm^-3','cm^-3', \
                                        'cm^-3','cm^-3','cm^-3', \
                                        'cm^-3','cm^-3', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1', \
                                        'J m^-2 s^-1 Hz^-1','']
            else:
              self.__nvar_atmo = 19
              self.__vars_atmo = ['ltau','T','Pg', \
                                  'Bx','By','Bz','vx', \
                                  'vy','vz','vturb','ne', \
                                  'nHT','nHa','nH(0)','nH(1)', \
                                  'nH(2)','nH(3)','nH(4)','np','f']
              self.__vars_atmo_units = ['','K','dyn cm^-2','G', \
                                        'G','G','km s^-1','km s^-1', \
                                        'km s^-1','km s^-1','cm^-3', \
                                        'cm^-3','cm^-3','cm^-3', \
                                        'cm^-3','cm^-3','cm^-3', \
                                        'cm^-3','cm^-3','']
            if debug: print(f'nvar_atmo {self.__nvar_atmo}')

            # Atmosphere pixel size                          v f_diff
            self.__s_atmo_c = self.__nz*self.__nvar_atmo*4 + 4
            if debug: print(f's_atmo_c {self.__s_atmo_c}')

            # Atmosphere size
            __s_atmo = self.__nx*self.__ny*self.__s_atmo_c
            if debug: print(f's_atmo {__s_atmo}')

            # Skip atmosphere
            f.seek(__s_atmo,1)

            # Get header for result
            self.__nvar = struct.unpack('i',f.read(4))[0]
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__jump_to_lambda = self.__jump_to_atmo + __s_atmo + 8
            f.seek(self.__nl*8,1)
            if debug: print(f'nvar {self.__nvar}')
            if debug: print(f'nl {self.__nl}')

            # Prepare lists with flags and nodes for inversion
            # variables
            self.__inv_flag = []
            self.__inv_node = []

            # Primordial variable list
            self.__vars = [self.__bx,self.__by,self.__bz, \
                           'f','T',self.__vx,self.__vy,self.__vz, \
                           'vturb','Pg','Re{J21}','Im{J21}', \
                           'Re{J22}','Im{J22}']
            self.__vars_units = [self.__bxu,self.__byu,self.__bzu, \
                                 '','K',self.__vxu,self.__vyu, \
                                 self.__vzu,'km s^-1','dyn cm^-2', \
                                 'J m^-2 s^-1 Hz^-1', \
                                 'J m^-2 s^-1 Hz^-1', \
                                 'J m^-2 s^-1 Hz^-1', \
                                 'J m^-2 s^-1 Hz^-1']

            # Additional variables in pixel result
            if self.__polarization:
                self.__vars_add = ['chi2_0','chi2', \
                                   'Io','Qo','Uo','Vo', \
                                   'If','Qf','Uf','Vf']
                self.__vars_add_units = ['','', \
                                         'J m^-2 s^-1 sr Hz^-1', \
                                         'J m^-2 s^-1 sr Hz^-1',
                                         'J m^-2 s^-1 sr Hz^-1',
                                         'J m^-2 s^-1 sr Hz^-1',
                                         'J m^-2 s^-1 sr Hz^-1',
                                         'J m^-2 s^-1 sr Hz^-1',
                                         'J m^-2 s^-1 sr Hz^-1',
                                         'J m^-2 s^-1 sr Hz^-1']
                self.__vars_stk = ['Io','Qo','Uo','Vo', \
                                   'If','Qf','Uf','Vf']

            else:
                self.__vars_add = ['chi2_0','chi2','Io','If']
                self.__vars_add_units = ['','', \
                                         'J m^-2 s^-1 sr Hz^-1', \
                                         'J m^-2 s^-1 sr Hz^-1']
                self.__vars_stk = ['Io','If']

            # Dictionary of var order
            self.__vars_dic = {0: self.__bx, \
                               1: self.__by, \
                               2: self.__bz, \
                               3: 'f', \
                               4: 'T', \
                               5: self.__vx, \
                               6: self.__vy, \
                               7: self.__vz, \
                               8: 'vturb', \
                               9: 'Pg', \
                              10: 'Re{J21}', \
                              11: 'Im{J21}', \
                              12: 'Re{J22}', \
                              13: 'Im{J22}'}

            # For each variable
            for ivar in range(self.__nvar):

                # Read if inverting
                self.__inv_flag.append( \
                        struct.unpack('i', f.read(4))[0] > 0)

                # Number of nodes
                self.__inv_node.append( \
                        struct.unpack('i', f.read(4))[0])

                if (debug):
                  print(f'ivar {ivar} flag {self.__inv_flag[-1]} ' + \
                        f'node {self.__inv_node[-1]}')

            # Remove from variable list the ones without nodes
            for ivar in range(self.__nvar-1,-1,-1):
                if self.__inv_node[ivar] < 1:
                    self.__vars.pop(ivar)
                    self.__vars_units.pop(ivar)

            #
            # Prepare jumps results related

            # Size of header
            __s_res_h = 4 + 4 + self.__nvar*8 + self.__nl*8

            # Size of results, Stokes part
            if self.__polarization:

                self.__s_res_c = 8 + 32*self.__nl

            else:

                self.__s_res_c = 8 + 8*self.__nl

            # For each variable, add nodes info
            for flag,node in zip(self.__inv_flag,self.__inv_node):

                # Nodes
                if node <= 0: continue

                # Inverting
                if flag:
                    self.__s_res_c += 12*node
                else:
                    self.__s_res_c += 8*node

            # Results size
            __s_res = self.__nx*self.__ny*self.__s_res_c

            # Jump to beginning of results
            self.__jump_to_res = self.__s_head + __s_atmo + __s_res_h

            # Try reading RF
            try_rf = True

        except struct.error:

            try_rf = False
            raise

        # Initialize sizes
        __s_rf_head = 0
        self.__s_rf_c = 0

        # If tryin reading response function
        if try_rf:

            # Try reading response function header
            try:

                # Skip results
               #f.seek(__s_res,1)
                f.seek(self.__jump_to_res + __s_res,0)

                # Read variables with RF
                self.__nvar_rf = struct.unpack('i',f.read(4))[0]
                if debug: print('nvar_rf',self.__nvar_rf)
                __s_rf_head = 4

                # Prepare number of varying nodes
                self.__inv_vnode = []

                # For each variable
                for ivar,flag in enumerate(self.__inv_flag):

                    if debug: print('\nGet var',ivar,flag)

                    # Skip not inverting
                    if not flag:
                        self.__inv_vnode.append(0)
                        if debug: print('No ivert')
                        continue

                    # Read index
                    jvar = struct.unpack('i',f.read(4))[0]
                    if debug: print('Get var index',jvar)

                    # Read varying nodes
                    self.__inv_vnode.append( \
                            struct.unpack('i',f.read(4))[0])
                    if debug: print('Get nodes',self.__inv_vnode[-1])

                    # Add sizes
                    __s_rf_head += 8

                    # If varying nodes
                    if self.__inv_vnode[-1] > 0:
                        if self.__polarization:
                            self.__s_rf_c += \
                               (4 + 16*self.__nl)*self.__inv_vnode[-1]
                        else:
                            self.__s_rf_c += \
                               (4 + 4*self.__nl)*self.__inv_vnode[-1]

                # Jump to RF
                self.__jump_to_rf = self.__jump_to_res + \
                                    __s_res + __s_rf_head
                self.__is_RF = True

            except struct.error:

                self.__is_RF = False
                self.__jump_to_rf = 0

            except:

                raise

        else:

            self.__jump_to_rf = 0
            self.__is_RF = False

        # Close file
        f.close()

        # Get real size
        real_size = os.path.getsize(self.__filename)
        if self.__is_RF:
            expectedsize = self.__jump_to_rf + \
                           self.__nx*self.__ny*self.__s_rf_c
        else:
            expectedsize = self.__jump_to_res + __s_res

        if real_size == expectedsize:
            self.__complete = True
        else:
            self.__complete = False
            msg = f'This is an incomplete file. \n' + \
                  f'Expected size {expectedsize}, ' + \
                  f'but got {real_size} instead'
            _error(msg,0)

        if debug:
            print('\nSizes:')
            print('Head',self.__s_head)
            print('J2 atmo',self.__jump_to_atmo)
            print('atmo c',self.__s_atmo_c)
            print('atmo',__s_atmo)
            print('J2 lamb',self.__jump_to_lambda)
            print('Head res',__s_res_h)
            print('Res c',self.__s_res_c)
            print('Res',__s_res)
            print('J2 res',self.__jump_to_res)
            print('Try RF?',try_rf)
            print('RF?',self.__is_RF)
            print('RF head',__s_rf_head)
            print('RF c',self.__s_rf_c)
            print('RF',self.__nx*self.__ny*self.__s_rf_c)
            print('J2 RF',self.__jump_to_rf)
            print('is RF',self.__is_RF)

            print('\n')

        #
        # Create memmaps

        # Wavelength always available
        self.__lam = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))

        # Only if the file is complete
        if self.__complete:

            #
            # Stokes block

            # Chi is always
            fields = [('chi',np.float32,(2))]

            # If polarized
            if self.__polarization:
                fields.append(('stokeso',np.float32,(4,self.__nl)))
                fields.append(('stokesf',np.float32,(4,self.__nl)))
            else:
                fields.append(('stokeso',np.float32,(1,self.__nl)))
                fields.append(('stokesf',np.float32,(1,self.__nl)))

            #
            # Nodes block
            for ipar,node in enumerate(self.__inv_node):

                # No nodes, skip
                if node <= 0: continue

                # Variable
                var = self.__vars_dic[ipar]

                # Add
                if self.__inv_flag[ipar]:
                    fields.append((var,np.float32,(3,node)))
                else:
                    fields.append((var,np.float32,(2,node)))

            # Data type
            self.__fit_dtype = np.dtype(fields)

            # Get memmap
            self.__fit = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__jump_to_res, \
                                   dtype=self.__fit_dtype, \
                                   shape=(self.__nx,self.__ny))

            #
            # Atmosphere block

            # Chi is always
            fields = [('atmos',np.float32,(self.__nvar_atmo, \
                                           self.__nz)), \
                      ('f',np.float32,(1))]

            # Data type
            self.__atmo_dtype = np.dtype(fields)

            # Get memmap
            self.__atmo = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__jump_to_atmo, \
                                    dtype=self.__atmo_dtype, \
                                    shape=(self.__nx,self.__ny))
            #
            # RF block
            if self.__is_RF:

                # Initialize
                fields = []

                # For each variable
                for jvar,flag in enumerate(self.__inv_flag):

                    # Skip not inverting
                    if not flag: continue

                    # Varying nodes
                    vnode = self.__inv_vnode[jvar]
                    if vnode < 1: continue

                    # Read heights
                    fields.append((f'H_{self.__vars_dic[jvar]}', \
                                   np.float32,(vnode)))

                # For each variable
                for jvar,flag in enumerate(self.__inv_flag):

                    # Skip not inverting
                    if not flag: continue

                    # Varying nodes
                    vnode = self.__inv_vnode[jvar]
                    nl = self.__nl

                    # No nodes
                    if vnode < 1: continue

                    # Read RF
                    vl = f'RF_{self.__vars_dic[jvar]}'
                    if self.__polarization:
                        fields.append((vl,np.float32, \
                                       (vnode,4,self.__nl)))
                    else:
                        fields.append((vl,np.float32, \
                                       (vnode,1,self.__nl)))

                # Data type
                self.__rf_dtype = np.dtype(fields)

                # Get memmap
                self.__rf = np.memmap(self.__filename, \
                                      mode='r', \
                                      offset=self.__jump_to_rf, \
                                      dtype=self.__rf_dtype, \
                                      shape=(self.__nx,self.__ny))

        # Return valid
        return True


    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_polarized(self):
        ''' Get if the inversion has polarization
        '''
        return self.__polarization

    def _get_vtype(self):
        ''' Get the reference frame of the velocity vector
        '''
        if self.__vtype == 0:
            return 'vertical'
        else:
            return 'LOS'

    def _get_btype(self):
        ''' Get the reference frame of the magnetic field vector
        '''
        if self.__btype == 0:
            return 'vertical'
        else:
            return 'LOS'

    def _get_jkqin(self):
        ''' Get if there are jkq tensors in the inversion
        '''
        return self.__jkqa

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nxyz(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self.__nx, self.__ny, self.__nz

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_nvar_atmo(self):
        ''' Get number of variables in model atmosphere
        '''
        return self.__nvar_atmo

    def _get_nvar(self):
        ''' Get number of variables in node data
        '''
        return self.__nvar

    def _get_dims(self):
        ''' Get number of positions in x, y, height, and
            wavelength axes
        '''
        return self.__nx, self.__ny, self.__nz, self.__nl

    def _get_vars(self):
        ''' Get variables with node data
        '''
        return self.__vars + self.__vars_add

    def _get_vars_fix(self):
        ''' Get variables without nodes
        '''
        return self.__vars_add

    def _get_vars_units(self):
        ''' Get variables with node data with units
        '''
        out = []
        for var,uni in zip(self.__vars,self.__vars_units):
            out.append(var+' ['+uni+']')
        for var,uni in zip(self.__vars_add,self.__vars_add_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_vars_fix_units(self):
        ''' Get variables without node data with units
        '''
        out = []
        for var,uni in zip(self.__vars_add,self.__vars_add_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_vars_atmo(self):
        ''' Get variables available in the atmospheric model
        '''
        return self.__vars_atmo

    def _get_vars_atmo_units(self):
        ''' Get variables available in the atmospheric model
            with units
        '''
        out = []
        for var,uni in zip(self.__vars_atmo,self.__vars_atmo_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = self.__lam.copy()
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def _get_column(self,ix,iy,minl=None,maxl=None,var=None):
        ''' Get result of inversion for a given pixel
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if var is None:
            ivar = self.__vars + self.__vars_add
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__vars + self.__vars_add:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars',1)
                    return None

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self.__lam.copy()

        # If the file is complete
        if self.__complete:

            # Try getting data
            try:

                # Chi_ori
                chi0 = self.__fit['chi'][ix,iy].copy()
                chi = chi0[1]
                chi0 = chi0[0]

                # Stokes
                stokes_ob = self.__fit['stokeso'][ix,iy].copy()
                stokes_fi = self.__fit['stokesf'][ix,iy].copy()

                # Nodes
                nodes = {}

                # For each variable
                for ipar,node in enumerate(self.__inv_node):

                    # No nodes
                    if node <= 0: continue

                    # Current variable
                    var = self.__vars_dic[ipar]

                    # Prepare space for this variable
                    nodes[var] = np.empty((3,node))

                    # Read
                    if self.__inv_flag[ipar]:
                        nodes[var] = self.__fit[var][ix,iy].copy()
                    else:
                        nodes[var][:2,:] = \
                                     self.__fit[var][ix,iy].copy()
                        nodes[var][2,:] = np.zeros((node))

            # Others
            except:
                raise

        # Incomplete file
        else:

            # Get column size
            bsiz = self.__s_res_c
            siz = bsiz//4

            # Need lambda?
            if minl is not None or maxl is not None:
                lam = self._get_lambda()

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__jump_to_res + iy*bsiz + \
                       self.__ny*ix*bsiz,0)

                # Chi_ori
                chi0 = struct.unpack('f',f.read(4))[0]
                chi = struct.unpack('f',f.read(4))[0]

                # Stokes
                if self.__polarization:
                    stokes_ob = np.array( \
                            struct.unpack('f'*self.__nl*4, \
                                          f.read(16*self.__nl))). \
                                        reshape((4,self.__nl))
                    stokes_fi = np.array( \
                            struct.unpack('f'*self.__nl*4, \
                                          f.read(16*self.__nl))). \
                                        reshape((4,self.__nl))
                else:
                    stokes_ob = np.array( \
                            struct.unpack('f'*self.__nl, \
                                          f.read(4*self.__nl))). \
                                        reshape((1,self.__nl))
                    stokes_fi = np.array( \
                            struct.unpack('f'*self.__nl, \
                                          f.read(4*self.__nl))). \
                                        reshape((1,self.__nl))

                # Nodes
                nodes = {}

                # For each variable
                for ipar,node in enumerate(self.__inv_node):

                    # No nodes
                    if node <= 0: continue

                    # Current variable
                    var = self.__vars_dic[ipar]

                    # Prepare space for this variable
                    nodes[var] = np.empty((3,node))

                    # Read H and value
                    nodes[var][0:2,:] = np.array( \
                            struct.unpack('f'*node*2, \
                                          f.read(8*node))). \
                                                reshape((2,node))
                    # Error?
                    if self.__inv_flag[ipar]:
                        nodes[var][2,:] = np.array( \
                                struct.unpack('f'*node, \
                                              f.read(4*node)))
                    else:
                        nodes[var][2,:] = 0e0

                # Close
                f.close()

            # Failed
            except struct.error:

                # If the file is complete, the error is more
                # severe, let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Close and warn
                    f.close()
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    chi0 = 0.
                    chi = 0.
                    if self.__polarization:
                        stokes_ob = np.zeros((4,self.__nl))
                        stokes_fi = np.zeros((4,self.__nl))
                    else:
                        stokes_ob = np.zeros((1,self.__nl))
                        stokes_fi = np.zeros((1,self.__nl))

                    # Nodes
                    nodes = {}

                    # For each variable
                    for ipar,node in enumerate(self.__inv_node):

                        # No nodes
                        if node <= 0: continue

                        # Current variable
                        var = self.__vars_dic[ipar]

                        # Generate zeros for this variable
                        nodes[var] = np.zeros((3,node))

            # Others
            except:
                raise

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            stokes_ob = stokes_ob[:,i:]
            stokes_fi = stokes_fi[:,i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            stokes_ob = stokes_ob[:,:i+1]
            stokes_fi = stokes_fi[:,:i+1]

        # Return column
        out = {}
        if 'chi2_0' in ivar:
            out['chi2_0'] = chi0
        if 'chi2' in ivar:
            out['chi2'] = chi
        if 'Io' in ivar:
            out['Io'] = stokes_ob[0,:]
        if 'Qo' in ivar:
            out['Qo'] = stokes_ob[1,:]
        if 'Uo' in ivar:
            out['Uo'] = stokes_ob[2,:]
        if 'Vo' in ivar:
            out['Vo'] = stokes_ob[3,:]
        if 'If' in ivar:
            out['If'] = stokes_fi[0,:]
        if 'Qf' in ivar:
            out['Qf'] = stokes_fi[1,:]
        if 'Uf' in ivar:
            out['Uf'] = stokes_fi[2,:]
        if 'Vf' in ivar:
            out['Vf'] = stokes_fi[3,:]
        if 'T' in ivar:
            out['T'] = nodes['T']
        if 'Pg' in ivar:
            out['Pg'] = nodes['Pg']
        if self.__bx in ivar:
            out[self.__bx] = nodes[self.__bx]
        if self.__by in ivar:
            out[self.__by] = nodes[self.__by]
        if self.__bz in ivar:
            out[self.__bz] = nodes[self.__bz]
        if self.__vx in ivar:
            out[self.__vx] = nodes[self.__vx]
        if self.__vy in ivar:
            out[self.__vy] = nodes[self.__vy]
        if self.__vz in ivar:
            out[self.__vz] = nodes[self.__vz]
        if 'vturb' in ivar:
            out['vturb'] = nodes['vturb']
        if self.__jkqa:
            if 'Re{J21}' in ivar:
                out['Re{J21}'] = nodes['Re{J21}']
            if 'Im{J21}' in ivar:
                out['Im{J21}'] = nodes['Im{J21}']
            if 'Re{J22}' in ivar:
                out['Re{J22}'] = nodes['Re{J22}']
            if 'Im{J22}' in ivar:
                out['Im{J22}'] = nodes['Im{J22}']

        return out

    def _get_column_atmo(self,ix,iy,minh=None,maxh=None,var=None):
        ''' Get model atmosphere for a given pixel
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if var is None:
            ivar = self.__vars_atmo
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__vars_atmo:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_atmo',1)
                    return None

        # If file is complete
        if self.__complete:

            # Try getting data
            try:

                # Read data
                col = self.__atmo['atmos'][ix,iy].copy()

                # Read f
                f_diff = self.__atmo['f'][ix,iy].copy()[0]

                # Check upper tau
                if col[0,0] <= 0.: col[0,0] = col[0,1]*1e-3

                # Logarithm
                col[0,:] = np.log10(col[0,:])

                # Adjust height
                if minh is not None:
                    i = np.argmin(np.absolute(col[0,:] - minh))
                    col = col[:,:i+1]
                if maxh is not None:
                    i = np.argmin(np.absolute(col[0,:] - maxh))
                    col = col[:,i:]

            except:
                raise

        # If incomplete file
        else:

            # Get column size
            bsiz = self.__s_atmo_c
            siz = bsiz//4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__jump_to_atmo + iy*bsiz + \
                       self.__ny*ix*bsiz,0)

                # Remove diffuse light factor from size
                bsiz -= 4
                siz -= 1

                # Read data
                col = np.array(struct.unpack('f'*siz, \
                                             f.read(bsiz))). \
                               reshape((self.__nvar_atmo,self.__nz))
                # Read f
                f_diff = struct.unpack('f',f.read(4))[0]

                # Close
                f.close()

                # Check upper tau
                if col[0,0] <= 0.: col[0,0] = col[0,1]*1e-3

                # Logarithm
                col[0,:] = np.log10(col[0,:])

                # Adjust height
                if minh is not None:
                    i = np.argmin(np.absolute(col[0,:] - minh))
                    col = col[:,:i+1]
                if maxh is not None:
                    i = np.argmin(np.absolute(col[0,:] - maxh))
                    col = col[:,i:]

            except struct.error:

                # If the file is complete, the error is more
                # severe, let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Close and warn
                    f.close()
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    col = np.zeros((self.__nvar_atmo,self.__nz))
                    f_diff = 0

            except:
                raise

        # Return column
        out = {}
        if 'ltau' in ivar:
            out['ltau'] = col[0,:]
        if 'T' in ivar:
            out['T'] = col[1,:]
        if 'Pg' in ivar:
            out['Pg'] = col[2,:]
        if 'Bx' in ivar:
            out['Bx'] = col[3,:]
        if 'By' in ivar:
            out['By'] = col[4,:]
        if 'Bz' in ivar:
            out['Bz'] = col[5,:]
        if 'vx' in ivar:
            out['vx'] = col[6,:]
        if 'vy' in ivar:
            out['vy'] = col[7,:]
        if 'vz' in ivar:
            out['vz'] = col[8,:]
        if 'vturb' in ivar:
            out['vturb'] = col[9,:]
        if 'ne' in ivar:
            out['ne'] = col[10,:]
        if 'nHT' in ivar:
            out['nHT'] = col[11,:]
        if 'nHa' in ivar:
            out['nHa'] = col[12,:]
        if 'nH(0)' in ivar:
            out['nH(0)'] = col[13,:]
        if 'nH(1)' in ivar:
            out['nH(1)'] = col[14,:]
        if 'nH(2)' in ivar:
            out['nH(2)'] = col[15,:]
        if 'nH(3)' in ivar:
            out['nH(3)'] = col[16,:]
        if 'nH(4)' in ivar:
            out['nH(4)'] = col[17,:]
        if 'np' in ivar:
            out['np'] = col[18,:]
        if 'f' in ivar:
            out['f'] = f_diff
        if self.__jkqa:
            if 'J10' in ivar:
                out['J10'] = col[19,:]
            if 'Re{J11}' in ivar:
                out['Re{J11}'] = col[20,:]
            if 'Im{J11}' in ivar:
                out['Im{J11}'] = col[21,:]
            if 'J20' in ivar:
                out['J20'] = col[22,:]
            if 'Re{J21}' in ivar:
                out['Re{J21}'] = col[23,:]
            if 'Im{J21}' in ivar:
                out['Im{J21}'] = col[24,:]
            if 'Re{J22}' in ivar:
                out['Re{J22}'] = col[25,:]
            if 'Im{J22}' in ivar:
                out['Im{J22}'] = col[26,:]

        return out

    def _get_column_rf(self,ix,iy,minl=None,maxl=None,var=None):
        ''' Get the response function for a given pixel
        '''

        # Is there RF?
        if not self.__is_RF:
           _error('This file does not have response functions',1)
           return None

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if var is None:
            ivar = self.__vars
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__vars:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars',1)
                    return None

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self.__lam.copy()

        # If complete file
        if self.__complete:

            # Output
            out = {}

            # Try getting data
            try:

                # Create dictionary
                red = {}

                # For each variable
                for jvar,flag in enumerate(self.__inv_flag):

                    # Skip not inverting
                    if not flag: continue

                    # Varying nodes
                    vnode = self.__inv_vnode[jvar]
                    if vnode < 1: continue

                    # labels
                    hl = f'H_{self.__vars_dic[jvar]}'
                    vl = f'RF_{self.__vars_dic[jvar]}'

                    # Read heights
                    red[jvar] = {'H': self.__rf[hl][ix,iy].copy(), \
                                 'RF': self.__rf[vl][ix,iy].copy()}

                    # Add to output?
                    if self.__vars_dic[jvar] in ivar:
                        out[self.__vars_dic[jvar]] = red[jvar]
            except:
                raise

        # Incomplete file
        else:

            # Get column size
            bsiz = self.__s_rf_c
            siz = bsiz//4

            # Output
            out = {}

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__jump_to_rf + iy*bsiz + \
                       self.__ny*ix*bsiz,0)

                # Create dictionary
                red = {}

                # For each variable
                for jvar,flag in enumerate(self.__inv_flag):

                    # Skip not inverting
                    if not flag: continue

                    # Varying nodes
                    vnode = self.__inv_vnode[jvar]
                    if vnode < 1: continue

                    # Read heights
                    red[jvar] = {'H': np.array( \
                                      struct.unpack('f'*vnode, \
                                                    f.read(4*vnode)))}

                # For each variable
                for jvar,flag in enumerate(self.__inv_flag):

                    # Skip not inverting
                    if not flag: continue

                    # Varying nodes
                    vnode = self.__inv_vnode[jvar]
                    nl = self.__nl

                    # No nodes
                    if vnode < 1: continue

                    # Read RF
                    if self.__polarization:
                        red[jvar]['RF'] = \
                            np.array( \
                            struct.unpack('f'*vnode*4*nl, \
                                          f.read(16*nl*vnode))). \
                                     reshape((vnode,4,self.__nl))
                    else:
                        red[jvar]['RF'] = \
                            np.array( \
                            struct.unpack('f'*vnode*nl, \
                                          f.read(4*nl*vnode))). \
                                     reshape((vnode,1,self.__nl))

                    # Add to output?
                    if self.__vars_dic[jvar] in ivar:
                        out[self.__vars_dic[jvar]] = red[jvar]

                # Close file
                f.close()

            # Failed
            except struct.error:

                # If the file is complete, the error is more
                # severe, let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Close and warn
                    f.close()
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Initialize
                    red = {}

                    # For each variable
                    for jvar,flag in enumerate(self.__inv_flag):

                        # Skip not inverting
                        if not flag: continue

                        # Varying nodes
                        vnode = self.__inv_vnode[jvar]
                        if vnode < 1: continue

                        # Read heights
                        red[jvar] = {'H': np.zeros(vnode)}

                    # Generate zeros
                    # For each variable
                    for jvar,flag in enumerate(self.__inv_flag):

                        # Skip not inverting
                        if not flag: continue

                        # Varying nodes
                        vnode = self.__inv_vnode[jvar]
                        nl = self.__nl

                        # Read RF
                        if self.__polarization:
                            red[jvar] = \
                                 {'RF': np.zeros((vnode,4,self.__nl))}
                        else:
                            red[jvar] = \
                                 {'RF': np.zeros((vnode,1,self.__nl))}

                        # Add to output?
                        if self.__vars_dic[jvar] in ivar:
                            out[ivar] = red[jvar]

            except:
                raise

        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            for var in out:
                out[var]['RF'] = out[var]['RF'][:,:,i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            for var in out:
                out[var]['RF'] = out[var]['RF'][:,:,:i+1]

        return out

    def _get_plane_chi(self):
        ''' Get full plane of chi2 and chi0
        '''

        # Complete file
        if self.__complete:

            # Fetch
            chi = self.__fit['chi'].copy()

        # Incomplete file
        else:

            # Get column size
            bsiz = self.__s_res_c

            # Before and after
            before = 0
            after = bsiz - 8

            # Initialize
            chi = np.empty((self.__nx,self.__ny,2))
            abort = False

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__jump_to_res,0)

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try to get chi
                    try:

                        chi[ix,iy,:] = np.array( \
                                        struct.unpack('ff',f.read(8)))

                    except struct.error:
                        if self.__complete:
                            raise
                        else:
                            msg = 'Could not read, may be due to ' + \
                                  'the file being not complete'
                            _error(msg,0)
                            chi[ix,iy:self.__ny,:] = 0.0
                            abort = True
                            break
                    except:
                        raise

                    # Jump to next
                    f.seek(after,1)

                # We are leaving
                if abort:
                    chi[ix+1:self.__nx,:,:] = 0.0
                    break

            # Close
            f.close()

        return chi


    def _get_plane_stk(self,il,var=None):
        ''' Get observed or fitted Stokes parameters for a given
            wavelength
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None


        # If var is not None
        if var is None:
            ivar = self.__vars_stk
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__vars_stk:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars',1)
                    return None

        # Variable names
        if self.__polarization:
            stko = ['Io','Qo','Uo','Vo']
            stkf = ['If','Qf','Uf','Vf']
        else:
            stko = ['Io']
            stkf = ['If']

        # Complete file
        if self.__complete:

            # Prepare output
            out = {}

            # Try getting data
            try:

                # For each stokes in observation
                for istk,stk in enumerate(stko):

                    # Output
                    if stk in ivar:

                        # Read
                        out[stk] = \
                                self.__fit['stokeso'][:,:,0,il].copy()

                # For each stokes in fit
                for istk,stk in enumerate(stkf):

                    # Output
                    if stk in ivar:

                        # Read
                        out[stk] = \
                                self.__fit['stokesf'][:,:,0,il].copy()
            except:
                raise

        # Incomplete file
        else:

            # Get column size
            bsiz = self.__s_res_c

            # Before and after
            before = 8
            left = il*4
            right = (self.__nl - il - 1)*4
            full = self.__nl*4
            if self.__polarization:
                after = bsiz - 8 - self.__nl*32
            else:
                after = bsiz - 8 - self.__nl*8

            # Prepare output
            out = {}
            for jvar in ivar:
                out[jvar] = np.empty((self.__nx,self.__ny))

            # Initialize
            abort = False

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__jump_to_res,0)

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Skip data before Stokes
                        f.seek(before,1)

                        # For each stokes in observation
                        for istk,stk in enumerate(stko):

                            # Output
                            if stk in out:

                                # Skip left
                                if left > 0: f.seek(left,1)

                                # Read
                                out[stk][ix,iy] = \
                                       struct.unpack('f',f.read(4))[0]
                                # Skip right
                                if right > 0: f.seek(right,1)

                            else:

                                # Skip
                                f.seek(full,1)

                        # For each stokes in observation
                        for istk,stk in enumerate(stkf):

                            # Output
                            if stk in out:

                                # Skip left
                                if left > 0: f.seek(left,1)

                                # Read
                                out[stk][ix,iy] = \
                                       struct.unpack('f',f.read(4))[0]
                                # Skip right
                                if right > 0: f.seek(right,1)

                            else:

                                # Skip
                                f.seek(full,1)

                    except struct.error:
                        if self.__complete:
                            raise
                        else:
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)
                            for stk in out:
                                out[stk][ix,iy:self.__ny] = 0.0
                            abort = True
                            break
                    except:
                        raise

                    # Jump to next
                    f.seek(after,1)

                # We are leaving
                if abort:
                    for stk in out:
                        out[stk][ix+1:self.__nx,:] = 0.0
                    break

            # Close
            f.close()

        return out


    def _get_plane_atmo(self,iz,var=None):
        ''' Get model atmosphere for a given optical depth
        '''

        # Valid?
        if not isinstance(iz, int) and not isinstance(iz, npint):
           _error('iz must be an integer',1)
           return None
        if iz < 0 or iz >= self.__nz:
           _error('The requested height index is out of bounds',1)
           return None


        # If var is not None
        if var is None:
            ivar = self.__vars_atmo
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__vars_atmo:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_atmo',1)
                    return None

        # Complete file
        if self.__complete:

            # Prepare output
            out = {}

            # Try getting data
            try:

                # For each variable in the model atmosphere
                for ibar,bar in enumerate(self.__vars_atmo):

                    # If diffuse factor
                    if bar == 'f' and bar in ivar:

                        # Read
                        out[bar] = self.__atmo['f'][:,:,0].copy()

                    # Any other
                    elif bar in ivar:

                        # Copy
                        out[bar] = \
                              self.__atmo['atmos'][:,:,ibar,iz].copy()

                        # If tau
                        if bar == 'ltau':
                            out[bar] = np.log10(out[bar])
            except:
                raise

        # Incomplete file
        else:

            # Get column size
            bsiz = self.__s_atmo_c

            # Before and after
            left = iz*4
            right = (self.__nz - iz - 1)*4
            full = self.__nz*4

            # Prepare output
            out = {}
            for jvar in ivar:
                out[jvar] = np.empty((self.__nx,self.__ny))

            # Initialize
            abort = False

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__jump_to_atmo,0)

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # For each variable in the model atmosphere
                        for ibar,bar in enumerate(self.__vars_atmo):

                            # If tau
                            if bar == 'ltau':

                                # If ltau in the list
                                if bar in out:

                                    # Left
                                    if left > 0: f.seek(left,1)

                                    # Get tau
                                    aux = \
                                       struct.unpack('f',f.read(4))[0]

                                    # Skip right
                                    if right > 0: f.seek(right,1)

                                    # If output tau
                                    if aux <= 0:
                                        out['ltau'] = -16.
                                    else:
                                        out['ltau'] = np.log10(aux)

                                # Not to output
                                else:

                                    # Skip
                                    f.seek(full,1)

                            # If diffuse factor
                            elif bar == 'f':

                                # If in output
                                if bar in out:

                                    # Read
                                    out[bar] = \
                                       struct.unpack('f',f.read(4))[0]

                                else:

                                    # Skip
                                    f.seek(4,1)

                            # Rest of variables
                            else:

                              # Output
                              if bar in out:

                                  # Skip left
                                  if left > 0: f.seek(left,1)

                                  # Read
                                  out[bar][ix,iy] = \
                                       struct.unpack('f',f.read(4))[0]

                                  # Skip right
                                  if right > 0: f.seek(right,1)

                              else:

                                  # Skip
                                  f.seek(full,1)

                    except struct.error:
                        if self.__complete:
                            raise
                        else:
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)
                            for bar in out:
                                out[bar][ix,iy:self.__ny] = 0.0
                            abort = True
                            break
                    except:
                        raise

                # We are leaving
                if abort:
                    for bar in out:
                        out[bar][ix+1:self.__nx,:] = 0.0
                    break

            # Close
            f.close()

        return out

    def _get_node(self,var):
        ''' Get cube result for a given variable
        '''

        # If var is not None
        if not isinstance(var, str):
            _error('The field var requires a string',1)
            return None
        if var not in self.__vars:
            _error('The requested variable ' + var + \
                   ' is not available, ' + \
                   'check with get_vars',1)
            return None

        # If the file is complete
        if self.__complete:

            # Index of this variable
            for key in self.__vars_dic:
                if var == self.__vars_dic[key]:
                    ipara = key
                    break

            # Check shape
            shape = self.__fit[var].shape

            # If we have errors
            if shape[2] == 3:

                # Just copy
                out = self.__fit[var].copy()

            # If not, we need to fake the errors
            else:

                # Initialize
                out = np.empty((self.__nx, \
                                self.__ny, \
                                3,shape[-1]))

                # Get known data
                out[:,:,:2,:] = self.__fit[var].copy()

                # Fill with zeros
                out[:,:,2,:] = 0e0

        # Incomplete file
        else:

            # Index of this variable
            for key in self.__vars_dic:
                if var == self.__vars_dic[key]:
                    ipara = key
                    break

            # Get column size
            bsiz = self.__s_res_c

            # Before and after
            if self.__polarization:
                before = 8 + self.__nl*32
            else:
                before = 8 + self.__nl*8

            # Initialize left and right
            left = 0.
            right = 0.

            # For each variable before
            for ipar in range(0,ipara):

                # No nodes
                if self.__inv_node[ipar] <= 0: continue

                # Inverting
                if self.__inv_flag[ipar]:
                    left += 12*self.__inv_node[ipar]
                else:
                    left += 8*self.__inv_node[ipar]

            # For each variable after
            for ipar in range(ipara+1,len(self.__inv_node)):

                # No nodes
                if self.__inv_node[ipar] <= 0: continue

                # Inverting
                if self.__inv_flag[ipar]:
                    right += 12*self.__inv_node[ipar]
                else:
                    right += 8*self.__inv_node[ipar]

            # Initialize output
            nnode = self.__inv_node[ipara]
            out = np.empty((self.__nx,self.__ny,3,nnode))

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__jump_to_res,0)

            # Get size of actual output
            if self.__inv_flag[ipara]:
                nsiz = 3
            else:
                nsiz = 2
                out[:,:,-1,:] = 0.

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Skip before and left
                        f.seek(before + left)

                        # Get nodes
                        for col in range(nsiz):
                            out[ix,iy,col,:] = np.array( \
                                    struct.unpack('f'*nnode, \
                                                  f.read(4*nnode)))
                    except struct.error:
                        if self.__complete:
                            raise
                        else:
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)
                            out[ix,iy:self.__ny,:,:] = 0.0
                            abort = True
                            break
                    except:
                        raise

                    # Jump to next
                    if right > 0: f.seek(right,1)

                # We are leaving
                if abort:
                    for stk in out:
                        out[ix+1:self.__nx,:,:,:] = 0.0
                    break

            # Close
            f.close()

        return out

    def _get_cube(self):
        ''' Get fit memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__fit

    def _get_cube_atmo(self):
        ''' Get atmosphere memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__atmo

######################################################################
######################################################################
######################################################################

class _cols_damp_15D():
    ''' Class to manage the 1.5D collisional rate files
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_na': \
          [None,'Get number of atoms in the original run'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height ' + \
                'dimensions'], \
         'get_dims': \
          [None,' Get number of positions in x, y, and height ' + \
                'axes and entries'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method. Note that ' + \
           'the raw data is not scaled and, if this a file with ' + \
           'collisions, you must multiply by 10^8 to get the ' + \
           'rate in s^-1']}
        # Collisions
        if self.__cols:
            self.__methods['get_type'] = \
              [None,'Get the type of collisional rates in the file']
            self.__methods['get_nentry'] = \
              [None,'Get number of collisional entries']
            self.__methods['get_column'] = \
              [{'ix': \
                 'Coordinate in the x dimension of the column ' + \
                 'to extract', \
                'iy': \
                 'Coordinate in the y dimension of the column ' + \
                 'to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the collisional rates for a particular ' + \
               'column']
            self.__methods['get_plane'] = \
              [{'iz': \
                 'Coordinate in the height dimension of the ' + \
                 'atmospheric parameters to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the collisional rates for a particular ' + \
               'height index for the whole field of view']
        # Damping
        elif self.__damp:
            self.__methods['get_nentry'] = \
              [None,'Get number of damping entries']
            self.__methods['get_column'] = \
              [{'ix': \
                 'Coordinate in the x dimension of the column ' + \
                 'to extract', \
                'iy': \
                 'Coordinate in the y dimension of the column ' + \
                 'to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the dampings for a particular column']
            self.__methods['get_plane'] = \
              [{'iz': \
                 'Coordinate in the height dimension of the ' + \
                 'atmospheric parameters to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the dampings for a particular height ' + \
               'index for the whole field of view']
        # Elastic rates
        else:
            self.__methods['get_nentry'] = \
              [None,'Get number of elastic rates entries']
            self.__methods['get_column'] = \
              [{'ix': \
                 'Coordinate in the x dimension of the column ' + \
                 'to extract', \
                'iy': \
                 'Coordinate in the y dimension of the column ' + \
                 'to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the elastic rates for a particular column']
            self.__methods['get_plane'] = \
              [{'iz': \
                 'Coordinate in the height dimension of the ' + \
                 'atmospheric parameters to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the elastic rates for a particular ' + \
               'height index for ' + \
               'the whole field of view']

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert collisional/damping file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')

            # Skip 3 first letter
            f.seek(3,0)

            # Read the fourth for the type of collisions
            self.__type = f.read(1).decode('utf-8')
            self.__cols = self.__type == 'l' or self.__type == 't'
            if not self.__cols:
                self.__damp = self.__type == 'a'

            # Get dimensions
            self.__na = int(struct.unpack('i',f.read(4))[0])
            self.__nx = int(struct.unpack('i',f.read(4))[0])
            self.__ny = int(struct.unpack('i',f.read(4))[0])
            self.__nz = int(struct.unpack('i',f.read(4))[0])
            self.__nn = int(struct.unpack('i',f.read(4))[0])

            # Close file
            f.close()

            # Size of head
            self.__head = 4*6

            # Size of column
            self.__column = 4*self.__nz*self.__nn

        except struct.error:
            return False
            raise
        except:
            raise

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size
        expectedsize = self.__head + self.__column*self.__nx*self.__ny

        # If complete
        if real_size == expectedsize:

            self.__complete = True

        # Incomplete file
        else:

            # Incomplete file
            self.__complete = False
            msg = f'This is an incomplete file. \n' + \
                  f'Expected size {expectedsize}, ' + \
                  f'but got {real_size} instead'
            _error(msg,0)

        #
        # Create memmaps

        # If complete
        if self.__complete:
            self.__data = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__head, \
                                    dtype=np.float32, \
                                    shape=(self.__nx, \
                                           self.__ny, \
                                           self.__nn, \
                                           self.__nz))

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_type(self):
        ''' Get type of collisions in this file
        '''
        if self.__type == 'l':
            _error('Collisions between level',0)
        elif self.__type == 't':
            _error('Collisions between terms',0)
        else:
            _error('Type not recognized',1)

    def _get_na(self):
        ''' Get number of atoms in the original run
        '''
        return self.__na

    def _get_nentry(self):
        ''' Get number of entries in file
        '''
        return self.__nn

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nxyz(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self.__nx, self.__ny, self.__nz

    def _get_dims(self):
        ''' Get number of positions in x, y, and height axes
            and entries
        '''
        return self.__nx, self.__ny, self.__nz, self.__nn

    def _get_column(self,ix,iy,ie=None):
        ''' Get collisional rates/dampings for a given pixel
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if ie is None:
            trim = False
        else:
            trim = True
            if not isinstance(ie, list):
                ivar = [ie]
            else:
                ivar = ie.copy()
            for evar in ivar:
                if not isinstance(evar, int):
                    _error('The field var requires integers',1)
                    return None
                if evar < 0 or evar >= self.__nn:
                    _error('The requested entry ' + evar + \
                           ' is out of limits, ' + \
                           'check with get_nentry',1)
                    return None

        # If complete file
        if self.__complete:

            # Try getting data
            try:

                # Fetch data
                col = self.__data[ix,iy,ivar,:].copy()

            except:
                raise

        # Incomplete file
        else:

            # Size column
            size = self.__column//4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*self.__column + \
                           self.__ny*ix*self.__column,0)

                # Read data
                col = np.array( \
                        struct.unpack('f'*size, \
                                      f.read(self.__column))). \
                               reshape((self.__nn,self.__nz))

                # Close
                f.close()

            except struct.error:

                # If the file is complete, the error is more severe,
                # let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Warn
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    for j in indx:
                        col = np.zeros((self.__nn,self.__nz))

            except:
                raise

            # If trimming
            if trim:
                todel = []
                for i in range(self.__nn-1,-1,-1):
                    if i not in ivar:
                        todel.append(i)
                col = np.delete(col,np.array(todel,dtype='int32'), \
                                axis=0)

        # If collisions, units factor
        if self.__cols: col *= 1e8

        # Return column
        return col

    def _get_plane(self,iz,ie=None):
        ''' Get collisional rates/dampings for a given height index
        '''

        # Valid?
        if not isinstance(iz, int) and not isinstance(iz, npint):
           _error('iz must be an integer',1)
           return None
        if iz < 0 or iz >= self.__nz:
           _error('The requested height index is out of bounds',1)
           return None

        # If var is not None
        if ie is None:
            trim = False
        else:
            trim = True
            if not isinstance(ie, list):
                ivar = [ie]
            else:
                ivar = ie.copy()
            for evar in ivar:
                if not isinstance(evar, int):
                    _error('The field var requires integers',1)
                    return None
                if evar < 0 or evar >= self.__nn:
                    _error('The requested entry ' + evar + \
                           ' is out of limits, ' + \
                           'check with get_nentry',1)
                    return None

        # Complete file
        if self.__complete:

            # Try getting data
            try:

                # Fetch
                out = self.__data[:,:,ivar,iz].copy()

            except:
                raise

        # Incomplete file
        else:

            # Before and after
            left = iz*4
            right = (self.__nz - iz - 1)*4
            full = self.__nz*4
            abort = False

            # Prepare output
            out = np.empty((self.__nx,self.__ny,self.__nn))

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__head,0)

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # For each variable
                        for ibar in range(self.__nn):

                            # If not in output
                            if ibar not in ivar:
                                f.seek(full,1)
                                continue

                            # Skip left
                            if left > 0: f.seek(left,1)

                            # Read
                            out[ix,iy,ibar] = \
                                    struct.unpack('f',f.read(4))[0]

                            # Skip right
                            if right > 0: f.seek(right,1)

                    except struct.error:

                        # If the file is complete, the error is
                        # more severe, let it crash
                        if self.__complete:
                            raise

                        # Incomplete file, may be missing data
                        else:

                            # Warn
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)

                            # Generate zeros
                            out[ix,iy:self.__ny,:] = 0.0
                            abort = True
                            break
                    except:
                        raise

                # We are leaving
                if abort:
                    out[ix+1:self.__nx,:,:] = 0.0
                    break

            # Close
            f.close()

            # If trimming
            if trim:
                todel = []
                for i in range(self.__nn-1,-1,-1):
                    if i not in ivar:
                        todel.append(i)
                out = np.delete(out,np.array(todel,dtype='int32'), \
                                axis=2)

        # If collisions, units factor
        if self.__cols: out *= 1e8

        return out

    def _get_cube(self):
        ''' Get data memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__data[:,:,:,::-1]

######################################################################
######################################################################
######################################################################

class _back_15D():
    ''' Class to manage the 1.5D background files
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e3
        #  Transformation to CGS
       #self.__unit_trans = 1e0/299792458e2

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height ' + \
                'dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, height, and ' + \
                'wavelength dimensions'], \
         'get_vars': \
          [None,'Get list of continuum variables'], \
         'get_vars_alias': \
          [None,'Get list of continuum variables with ' + \
                'their aliases'], \
         'get_vars_units': \
          [None,'Get list of continuum variables with their ' + \
                'corresponding units'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_alias()}'}, \
           'Extract the continuum variables at a ' + \
           'particular column'], \
          'get_plane': \
          [{'iz': \
             'Coordinate in the z dimension of the plane ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output ' + \
             '(see the available ones with get_vars_alias()}'}, \
           'Extract the continuum variables at a particular plane'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method. Note that ' + \
           'the raw data is not scaled to the SI units and ' + \
          f'you need to multiply by {self.__unit_trans} to get ' + \
           'the correct emissivity. You can get this number ' + \
           'with the method get_scale()'], \
          'get_scale': \
          [None,f'Get the scale factor between the raw data and ' + \
           'SI units']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert background continuum file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)
            # Dimensions
            self.__nx = struct.unpack('i',f.read(4))[0]
            self.__ny = struct.unpack('i',f.read(4))[0]
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__nd = struct.unpack('i',f.read(4))[0]

        except struct.error:
            return False
            raise
        except:
            raise

        # Variables
        self.__nvar = 3
        self.__vars = [r'$\eta_{\rm c}$',r'$\sigma_{\rm c}$', \
                       r'$\epsilon{\rm c}$']
        self.__alias = ['eta','sig','eps']
        self.__vars_units = ['[m$^{-1}$]','[m$^{-1}$]', \
                          '[J m$^{-3}$ s$^{-1}$ sr$^{-1}$ Hz$^{-1}$]']

        # Sizes
        self.__jump_to_lambda = 6*4
        self.__head = self.__jump_to_lambda + self.__nl*8
        self.__column = 4*self.__nz*self.__nl*self.__nvar

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size
        expectedsize = self.__head + self.__column*self.__nx*self.__ny

        # If complete
        if real_size == expectedsize:

            self.__complete = True

        # Incomplete file
        else:

            # Incomplete file
            self.__complete = False
            msg = f'This is an incomplete file. \n' + \
                  f'Expected size {expectedsize}, ' + \
                  f'but got {real_size} instead'
            _error(msg,0)

        # Get memmaps
        self.__omg = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        if self.__complete:
            self.__data = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__head, \
                                    dtype=np.float32, \
                                    shape=(self.__nx, \
                                           self.__ny, \
                                           self.__nz, \
                                           self.__nvar, \
                                           self.__nl))

        # Return valid
        return True


    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nxyz(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self.__nx, self.__ny, self.__nz

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in x, y, height, and
            wavelength axes
        '''
        return self.__nx, self.__ny, self.__nz, self.__nl

    def _get_vars(self):
        ''' Get variables with node data
        '''
        return self.__vars

    def _get_vars_alias(self):
        ''' Get variables and their alias
        '''
        out = []
        for var,alias in zip(self.__vars,self.__alias):
            out.append(var+' -> ',alias)
        return out

    def _get_vars_units(self):
        ''' Get variables with node data with units
        '''
        out = []
        for var,uni in zip(self.__vars,self.__vars_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_scale(self):
        ''' Return the multiplicative factor to get the
            emissivity in SI units
        '''
        return self.__unit_trans

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except struct.error:
            raise
        except:
            raise

    def _get_column(self,ix,iy,minl=None,maxl=None,var=None):
        ''' Get background quantities for a given pixel
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # Complete file
        if self.__complete:

            # Try getting data
            try:

                # Fetch data
                col = self.__data[ix,iy,:,:,::-1].copy()

            except:
                raise

        # Incomplete file
        else:

            # Size
            siz = self.__column//4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*self.__column + \
                           self.__ny*ix*self.__column,0)

                # Read data
                col = np.array( \
                        struct.unpack('f'*siz, \
                                      f.read(self.__column))). \
                               reshape((self.__nz,self.__nvar, \
                                        self.__nl))[:,:,::-1]

                # Close
                f.close()

            except struct.error:

                # If the file is complete, the error is more
                # severe, let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Close and warn
                    f.close()
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    col = np.zeros((self.__nz,self.__nvar,self.__nl))

            except:
                raise

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            col = col[:,:,i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            col = col[:,:,:i+1]

        # Return column
        out = {}
        if 'eta' in ivar:
            out['eta'] = col[:,0,:]
        if 'sig' in ivar:
            out['sig'] = col[:,1,:]
        if 'eps' in ivar:
            out['eps'] = col[:,2,:]*self.__unit_trans

        return out

    def _get_plane(self,iz,minl=None,maxl=None,var=None):
        ''' Get background quantities for a given plane
        '''

        # Valid?
        if not isinstance(iz, int) and not isinstance(iz, npint):
           _error('ix must be an integer',1)
           return None
        if iz < 0 or iz >= self.__nz:
           _error('The requested plane is out of bounds',1)
           return None

        # If var is not None
        if var is None:
            ivar = self.__alias
        else:
            if not isinstance(var, list):
                ivar = [var]
            else:
                ivar = var.copy()
            for evar in ivar:
                if not isinstance(evar, str):
                    _error('The field var requires strings',1)
                    return None
                if evar not in self.__alias:
                    _error('The requested variable ' + evar + \
                           ' is not available, ' + \
                           'check with get_vars_alias',1)
                    return None

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # Complete file
        if self.__complete:

            # Try getting data
            try:

                # Fetch
                col = self.__data[:,:,iz,:,::-1].copy()

            except:
                raise

        # Incomplete file
        else:

            # Size
            left = iz*self.__nvar*self.__nl*4
            right = (self.__nz - iz - 1)*self.__nvar*self.__nl*4
            siz = self.__nvar*self.__nl
            bsiz = siz*4
            abort = False

            # Prepare reading
            col = np.empty((self.__nx,self.__ny,self.__nvar, \
                            self.__nl))

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__head,0)

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Skip left
                        if left > 0: f.seek(left,1)

                        # Read
                        col[ix,iy,:,:] = np.array( \
                                struct.unpack('f'*siz, \
                                              f.read(bsiz))). \
                                       reshape((self.__nvar, \
                                                self.__nl))[:,::-1]

                        # Skip right
                        if right > 0: f.seek(right,1)

                    except struct.error:

                        # If the file is complete, the error is
                        # more severe, let it crash
                        if self.__complete:
                            raise

                        # Incomplete file, may be missing data
                        else:

                            # Warn
                            msg = 'Could not read, may be due ' + \
                                  'to the file being not complete'
                            _error(msg,0)

                            # Generate zeros
                            col[ix,iy:self.__ny,:,:] = 0.0
                            abort = True
                            break
                    except:
                        raise

                # We are leaving
                if abort:
                    col[ix+1:self.__nx,:,:,:] = 0.0
                    break

            # Close
            f.close()

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            col = col[:,:,:,i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            col = col[:,:,:,:i+1]

        # Return column
        out = {}
        if 'eta' in ivar:
            out['eta'] = col[:,:,0,:]
        if 'sig' in ivar:
            out['sig'] = col[:,:,1,:]
        if 'eps' in ivar:
            out['eps'] = col[:,:,2,:]*self.__unit_trans

        return out

    def _get_cube(self):
        ''' Get background memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__data[:,:,:,:,::-1]

######################################################################
######################################################################
######################################################################

class _pop_dep_15D():
    ''' Class to manage the 1.5D population/departure files
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
       #self.__unit_trans = 1e6
        #  Remain cgs
        self.__unit_trans = 1e0

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_nentry': \
          [None,'Get number of entries'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height ' + \
                'dimensions'], \
         'get_dims': \
          [None,' Get number of positions in x, y, and height ' + \
                'axes and entries'], \
         'get_column': \
          [{'ix': \
            'Coordinate in the x dimension of the column to ' + \
            'extract', \
            'iy': \
            'Coordinate in the y dimension of the column to ' + \
            'extract', \
            'ie': \
            'List of entries to include in the output'}, \
            'Extract the populations/departures for a ' + \
            'particular column'], \
          'get_plane': \
          [{'iz': \
            'Coordinate in the height dimension of the ' + \
            'atmospheric parameters to extract', \
            'ie': \
            'List of entries to include in the output'}, \
            'Extract the populations/departures for a ' + \
            'particular height index for ' + \
            'the whole field of view'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method.']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert populations/departure coeffs. file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')

            # Skip 3 first letter
            f.seek(3,0)

            # Read the fourth for the type of file
            self.__type = f.read(1).decode('utf-8')

            # Get dimensions
            self.__nx = int(struct.unpack('i',f.read(4))[0])
            self.__ny = int(struct.unpack('i',f.read(4))[0])
            self.__nz = int(struct.unpack('i',f.read(4))[0])
            self.__nn = int(struct.unpack('i',f.read(4))[0])

            # Close file
            f.close()

            # Size of head
            self.__head = 4*5

            # Size of column
            self.__column = 4*self.__nz*self.__nn


        except struct.error:
            return False
            raise
        except:
            raise

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size
        expectedsize = self.__head + self.__column*self.__nx*self.__ny

        # If complete
        if real_size == expectedsize:

            self.__complete = True

        # Incomplete file
        else:

            # Incomplete file
            self.__complete = False
            msg = f'This is an incomplete file. \n' + \
                  f'Expected size {expectedsize}, ' + \
                  f'but got {real_size} instead'
            _error(msg,0)

        # Create memmaps
        if self.__complete:
            self.__data = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__head, \
                                    dtype=np.float32, \
                                    shape=(self.__nx, \
                                           self.__ny, \
                                           self.__nz, \
                                           self.__nn))

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_type(self):
        ''' Get type of population/departure file
        '''
        if self.__type == 'p':
            _error('Number density',0)
        elif self.__type == 'b':
            _error('Departure coefficient',0)
        else:
            _error('Type not recognized',1)

    def _get_nentry(self):
        ''' Get number of entries in file
        '''
        return self.__nn

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nz(self):
        ''' Get number of heights
        '''
        return self.__nz

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nxyz(self):
        ''' Get number of positions in x, y, and height axes
        '''
        return self.__nx, self.__ny, self.__nz

    def _get_dims(self):
        ''' Get number of positions in x, y, and height axes
            and entries
        '''
        return self.__nx, self.__ny, self.__nz, self.__nn

    def _get_column(self,ix,iy,ie=None):
        ''' Get populations/departures for a given pixel
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        # If var is not None
        if ie is None:
            trim = False
        else:
            trim = True
            if not isinstance(ie, list):
                ivar = [ie]
            else:
                ivar = ie.copy()
            for evar in ivar:
                if not isinstance(evar, int):
                    _error('The field var requires integers',1)
                    return None
                if evar < 0 or evar >= self.__nn:
                    _error('The requested entry ' + evar + \
                           ' is out of limits, ' + \
                           'check with get_nentry',1)
                    return None

        # Complete file
        if self.__complete:

            # Try getting data
            try:

                # Read data
                col = self.__data[ix,iy,:,ivar].copy()

            except:
                raise

        # Incomplete file
        else:

            # Size column
            size = self.__column//4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*self.__column + \
                           self.__ny*ix*self.__column,0)

                # Read data
                col = np.array( \
                        struct.unpack('f'*size, \
                                      f.read(self.__column))). \
                               reshape((self.__nz,self.__nn))

                # Close
                f.close()

            except struct.error:

                # If the file is complete, the error is more severe,
                # let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Warn
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    for j in indx:
                        col = np.zeros((self.__nz,self.__nn))

            except:
                raise

            # If trimming
            if trim:
                todel = []
                for i in range(self.__nn-1,-1,-1):
                    if i not in ivar:
                        todel.append(i)
                col = np.delete(col,np.array(todel,dtype='int32'), \
                                axis=1)

        # Units
        if self.__type == 'p':
            col *= self.__unit_trans

        # Return column
        return col

    def _get_plane(self,iz,ie=None):
        ''' Get populations/departures for a given height index
        '''

        # Valid?
        if not isinstance(iz, int) and not isinstance(iz, npint):
           _error('iz must be an integer',1)
           return None
        if iz < 0 or iz >= self.__nz:
           _error('The requested height index is out of bounds',1)
           return None

        # If var is not None
        if ie is None:
            trim = False
        else:
            trim = True
            if not isinstance(ie, list):
                ivar = [ie]
            else:
                ivar = ie.copy()
            for evar in ivar:
                if not isinstance(evar, int):
                    _error('The field var requires integers',1)
                    return None
                if evar < 0 or evar >= self.__nn:
                    _error('The requested entry ' + evar + \
                           ' is out of limits, ' + \
                           'check with get_nentry',1)
                    return None

        # Complete file
        if self.__complete:

            # Try getting data
            try:

                # Read
                out = self.__data[:,:,iz,ivar].copy()

            except:
                raise

        # Incomplete file
        else:

            # Before and after
            left = iz*4*self.__nn
            right = (self.__nz - iz - 1)*4*self.__nn
            abort = False

            # Prepare output
            out = np.empty((self.__nx,self.__ny,self.__nn))

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this block
            f.seek(self.__head,0)

            # Run over columns
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Skip left
                        if left > 0: f.seek(left,1)

                        # Read
                        out[ix,iy,:] = np.array( \
                                struct.unpack('f'*self.__nn, \
                                              f.read(4*self.__nn)))
                        # Skip right
                        if right > 0: f.seek(right,1)

                    except struct.error:

                        # If the file is complete, the error is
                        # more severe, let it crash
                        if self.__complete:
                            raise

                        # Incomplete file, may be missing data
                        else:

                            # Warn
                            msg = 'Could not read, may be due to ' + \
                                  'the file being not complete'
                            _error(msg,0)

                            # Generate zeros
                            out[ix,iy:self.__ny,:] = 0.0
                            abort = True
                            break
                    except:
                        raise

                # We are leaving
                if abort:
                    out[ix+1:self.__nx,:,:] = 0.0
                    break

            # Close
            f.close()

            # If trimming
            if trim:
                todel = []
                for i in range(self.__nn-1,-1,-1):
                    if i not in ivar:
                        todel.append(i)
                out = np.delete(out,np.array(todel,dtype='int32'), \
                                axis=2)

        # Units
        if self.__type == 'p':
            out *= self.__unit_trans

        return out

    def _get_cube(self):
        ''' Get full memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__data

######################################################################
######################################################################
######################################################################

class _stokes_CLE():
    ''' Class to manage emergent Stokes parameters from CLE synthesis
    '''

    def __init__(self,filename):
        ''' Initialize class
        '''

        # Store filename
        self.__filename = filename

        # Get header
        if not self.__head(): return None

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        # Method
        self.__methods = {
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_type': \
          [None,'Get type of model'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, and ' + \
                'wavelength dimensions'], \
         'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
         'get_geometry': \
          [None,'Get geometric data'], \
         'get_stokesi_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract Stokes I at a particular column'], \
         'get_stokesq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the ' + \
             'column to extract', \
            'iy': \
             'Coordinate in the y dimension of the ' + \
             'column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular column'], \
         'get_stokesu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular column'], \
         'get_stokesv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular column'], \
         'get_linear_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a ' + \
           'particular column'], \
         'get_stokes_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column ' + \
             'to extract', \
            'iy': \
             'Coordinate in the y dimension of the column ' + \
             'to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular column'], \
          'get_plane_stk': \
          [{'il': \
             'Coordinate in the wavelength dimension of the ' + \
             'Stokes parameters to extract', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract Stokes parameters at a given wavelength ' + \
           'position for the whole field of view'], \
         'get_stokesi_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract Stokes I at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokesv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_linear_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a ' + \
           'particular wavelength index for the whole ' + \
           'field of view'], \
         'get_stokes_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular ' + \
           'wavelength index for the whole field of view'], \
          'get_cube': \
          [None,f'Get a memmap to the whole data. The file must ' + \
           'be complete to use this method. Note that ' + \
           'the raw data is not scaled to the SI units and ' + \
          f'you need to multiply by {self.__unit_trans}. ' + \
           'You can get this number with the method get_scale()'], \
          'get_scale': \
          [None,f'Get the scale factor between the raw data and ' + \
           'SI units'] \
           }

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D emergence file head
        '''
        try:

            # Get header data
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Mode
            self.__mode = struct.unpack('i',f.read(4))[0]

            # Wavelength
            self.__nl = struct.unpack('i',f.read(4))[0]

            # Lambda
            self.__jump_to_lambda = 4*3
            f.seek(self.__nl*8,1)

            #
            # Geometry
            #
            self.__jump_to_geometry = self.__jump_to_lambda + \
                                      8*self.__nl

            # Cartesian
            if self.__mode == 0:
                self.__nx = struct.unpack('i',f.read(4))[0]
                f.seek(self.__nx*8,1)
                self.__ny = struct.unpack('i',f.read(4))[0]
                f.seek(self.__ny*8,1)
                self.__head = self.__jump_to_geometry + \
                              8*(1 + self.__nx + self.__ny)
            # Slab or not cartesian
            else:
                self.__nx = struct.unpack('i',f.read(4))[0]
                self.__ny = 1
                self.__head = self.__jump_to_geometry + 4 + \
                              self.__nx*8*2

            # Close file
            f.close()

            # Column size
            self.__column_size = 16*self.__nl

        except struct.error:
            raise
        except:
            raise

        # Get real size
        real_size = os.path.getsize(self.__filename)

        # Expected size cartesian
        if self.__mode == 0:
            expectedsize = self.__head + \
                           self.__nl*self.__nx*self.__ny*4
        # Expected size slab or non-cartesian
        else:
            expectedsize = self.__head + \
                           self.__nx*self.__ny*(16 + \
                           self.__nl*4*8)

        # If intensity size
        if real_size == expectedsize:

            self.__complete = True

        # Incomplete file
        else:

            self.__complete = False
            msg = f'I have guessed that this is ' + \
                  f'an incomplete file'
            _error(msg,0)

        # Create memmaps
        self.__omg = np.memmap(self.__filename, \
                               mode='r', \
                               offset=self.__jump_to_lambda, \
                               dtype=np.float64, \
                               shape=(self.__nl))
        # Cartesian
        if self.__mode == 0:
            self.__geom_dtype = np.dtype((  \
                                      ('nx',np.int32,(1)), \
                                      ('x',np.float64,(self.__nx)), \
                                      ('ny',np.int32,(1)), \
                                      ('y',np.float64,(self.__ny))))
            self.__geom = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__jump_to_geometry, \
                                    dtype=self.__geom_dtype, \
                                    shape=(1))
        # Slab or non-cartesian
        else:
            self.__geom_dtype = np.dtype(( \
                                      ('n',np.int32,(1)), \
                                      ('xy',np.float64,(self.__nx,))))
            self.__geom = np.memmap(self.__filename, \
                                    mode='r', \
                                    offset=self.__jump_to_geometry, \
                                    dtype=self.__geom_dtype, \
                                    shape=(1))
        if self.__complete:
            self.__stk = np.memmap(self.__filename, \
                                   mode='r', \
                                   offset=self.__head, \
                                   dtype=np.float32, \
                                   shape=(self.__nx, \
                                          self.__ny, \
                                          4,self.__nl))

        # Return valid
        return True

    def _get_filename(self):
        ''' Get the name of the read file
        '''
        return self.__filename

    def _get_type(self):
        ''' Get the type of model
        '''
        if self.__mode == 0:
            return 'Cartesian'
        elif self.__mode == 1:
            return 'Slab'
        elif self.__mode == 2:
            return 'Non-Cartesian'
        else:
            return 'Unknown'

    def _get_polarized(self):
        ''' Get if there is polarization
        '''
        return self.__mode > 1

    def _get_nx(self):
        ''' Get number of positions in x axis
        '''
        return self.__nx

    def _get_ny(self):
        ''' Get number of positions in y axis
        '''
        return self.__ny

    def _get_nxy(self):
        ''' Get number of positions and x and y axes
        '''
        return self.__nx, self.__ny

    def _get_nl(self):
        ''' Get number of wavelengths
        '''
        return self.__nl

    def _get_dims(self):
        ''' Get number of positions in x, y, and wavelength axes
        '''
        return self.__nx, self.__ny, self.__nl

    def _get_scale(self):
        ''' Return the multiplicative factor to get SI units
        '''
        return self.__unit_trans

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            lam = 1e2/self.__omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            return lam
        except:
            raise

    def __get_geometry(self):
        ''' Return data about geometry
        '''
        try:
            # Cartesian
            if self.__mode == 0:
                x = self.__geom['x'][0,:].copy()
                y = self.__geom['y'][0,:].copy()
            # Slab or non-cartesian
            else:
                x = self.__geom['xy'][0,:,0].copy()
                y = self.__geom['xy'][0,:,1].copy()
            return x,y
        except:
            raise

    def __get_gen_column(self,ix,iy,minl=None,maxl=None, \
                         fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters column
        '''

        # Output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = 1e2/self.__omg[::-1]

        # Complete file
        if self.__complete:

            # Try getting data
            try:

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = self.__stk[ix,iy,0,::-1].copy()

                    # Out?
                    if 0 in indx:
                        out[0] = stkI

                # Q, U, and V
                for j in range(1,4):

                    # To output
                    if j in indx:

                        # Read Stokes
                        out[j] = self.__stk[ix,iy,j,::-1].copy()

                    # Manage units
                    if fractional and j in indx:
                        out[j] /= stkI

            # Fail
            except:
                raise

        # Incomplete file
        else:

            # Get size to read
            siz = self.__nl
            bsiz = siz*4

            # Try getting data
            try:

                # Open file
                f = open(self.__filename,'rb')

                # Seek first data points for this column
                f.seek(self.__head + iy*self.__column_size + \
                       self.__ny*ix*self.__column_size,0)

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = np.array(struct.unpack('f'*siz, \
                                                  f.read(bsiz)))[::-1]

                    # Out?
                    if 0 in indx:
                        out[0] = stkI

                # No intensity
                else:

                    # Skip
                    f.seek(bsiz,1)

                # Q, U, and V
                for j in range(1,4):

                    # To output
                    if j in indx:

                        # Read Stokes
                        out[j] = np.array( \
                                struct.unpack('f'*siz, \
                                              f.read(bsiz)))[::-1]
                    else:

                        # Skip
                        f.seek(bsiz,1)

                    # Manage units
                    if fractional and j in indx:
                        out[j] /= stkI

            # Failed
            except struct.error:

                # If the file is complete, the error is more severe,
                # let it crash
                if self.__complete:
                    raise

                # Incomplete file, may be missing data
                else:

                    # Warn
                    msg = 'Could not read, may be due to the ' + \
                          'file being not complete'
                    _error(msg,0)

                    # Generate zeros
                    for j in indx:
                        out[j] = np.zeros((self.__nl))

            # Others
            except:
                raise

            # Close file
            f.close()

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            for j in indx:
                out[j] = out[j][i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            for j in indx:
                out[j] = out[j][:i+1]

        # Units
        if fractional:
            if 0 in indx:
                out[0] *= self.__unit_trans
        else:
            for j in indx:
                out[j] *= self.__unit_trans

        # Return
        return out

    def _get_stokesi_column(self,ix,iy,minl=None,maxl=None):
        ''' Get intensity profile at a given column
        '''

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl,False,[0])[0]

    def _get_stokesq_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes Q profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[1])[1]

    def _get_stokesu_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes U profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix, int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy, int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[2])[2]

    def _get_stokesv_column(self,ix,iy,minl=None,maxl=None, \
                            fractional=False):
        ''' Get Stokes V profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix,int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy,int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        return self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[3])[3]

    def _get_linear_column(self,ix,iy,minl=None,maxl=None, \
                           fractional=False):
        ''' Get Stokes V profile at a given column
        '''

        # Mode?
        if self.__mode == 1:
           return np.zeros(self.__nl)

        # Valid?
        if not isinstance(ix,int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy,int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        qu = self.__get_gen_column(ix,iy,minl,maxl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_column(self,ix,iy,minl=None,maxl=None, \
                           fractional=False):
        ''' Get Stokes parameter at a given column
        '''

        # Mode?
        if self.__mode == 1:
           _error('The file is only intensity',1)
           return None

        # Valid?
        if not isinstance(ix,int) and not isinstance(ix, npint):
           _error('ix must be an integer',1)
           return None
        if not isinstance(iy,int) and not isinstance(iy, npint):
           _error('iy must be an integer',1)
           return None
        if ix < 0 or iy < 0 or ix >= self.__nx or iy >= self.__ny:
           _error('The requested column is out of bounds',1)
           return None

        iquv = self.__get_gen_column(ix,iy,minl,maxl, \
                                     fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))


    def __get_gen_plane(self,il,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters plane
        '''

        # Complete file
        if self.__complete:

            # Output
            out = [None,None,None,None]

            # Try getting data
            try:

                # Intensity
                if 0 in indx or fractional:

                    # Get intensity
                    stkI = self.__stk[:,:,0,il].copy()

                    # Out?
                    if 0 in indx:
                        out[0] = stkI

                # Q, U, and V
                for j in range(1,4):

                    # To output
                    if j in indx:

                        # Get Stokes
                        out[j] = self.__stk[:,:,j,il].copy()

                    # Manage units
                    if fractional and j in indx:
                        out[j][ix,iy] /= stkI
            except:
                raise

        # Incomplete file
        else:

            # Get size to read
            left = il*4
            right = (self.__nl - il - 1)*4
            full = self.__nl*4
            abort = False

            # Output
            out = [None,None,None,None]

            # For each index requested
            for j in indx:
                out[j] = np.empty((self.__nx,self.__ny))

            # Open file
            f = open(self.__filename,'rb')

            # Seek to data
            f.seek(self.__head,0)

            # For each column
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Try getting data
                    try:

                        # Intensity
                        if 0 in indx or fractional:

                            # Get intensity
                            if left > 0: f.seek(left,1)
                            stkI = struct.unpack('f',f.read(4))[0]
                            if right > 0: f.seek(right,1)

                            # Out?
                            if 0 in indx:
                                out[0][ix,iy] = stkI

                        # No intensity
                        else:

                            # Skip
                            f.seek(full,1)

                        # Q, U, and V
                        for j in range(1,4):

                            # To output
                            if j in indx:

                                # Get Stokes
                                if left > 0: f.seek(left,1)
                                out[j][ix,iy] = \
                                       struct.unpack('f',f.read(4))[0]
                                if right > 0: f.seek(right,1)

                            # No output
                            else:

                                # Skip
                                f.seek(full,1)

                            # Manage units
                            if fractional and j in indx:
                                out[j][ix,iy] /= stkI

                    # Reading error
                    except struct.error:

                        # If the file is complete, the error
                        # is more severe, let it crash
                        if self.__complete:
                            raise

                        # Incomplete file, may be missing data
                        else:

                            # Warn
                            msg = 'Could not read, may be ' + \
                                  'due to the file being not complete'
                            _error(msg,0)

                            # Generate zeros
                            for j in indx:
                                out[j][ix,iy:self.__ny] = 0.0
                            abort = True
                            break

                    except:
                        raise

                # We are leaving
                if abort:
                    for j in indx:
                        out[j][ix+1:self.__nx,:] = 0.0
                    break

            # Close file
            f.close()

        # Units
        if fractional:
            if 0 in indx:
                out[0] *= self.__unit_trans
        else:
            for j in indx:
                out[j] *= self.__unit_trans

        # Return
        return out

    def _get_stokesi_plane(self,il):
        ''' Get intensity profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,False,[0])[0]

    def _get_stokesq_plane(self,il,fractional=False):
        ''' Get Stokes Q profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,fractional,[1])[1]

    def _get_stokesu_plane(self,il,fractional=False):
        ''' Get Stokes U profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,fractional,[2])[2]

    def _get_stokesv_plane(self,il,fractional=False):
        ''' Get Stokes V profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        return self.__get_gen_plane(jl,fractional,[3])[3]

    def _get_linear_plane(self,il,fractional=False):
        ''' Get Stokes V profile at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        qu = self.__get_gen_plane(jl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_plane(self,il,fractional=False):
        ''' Get Stokes profiles at a given wavelength for the FoV
        '''

        # Valid?
        if not isinstance(il, int) and not isinstance(il, npint):
           _error('il must be an integer',1)
           return None
        if il < 0 or il >= self.__nl:
           _error('The requested wavelength is out of bounds',1)
           return None

        jl = self.__nl - il - 1
        iquv = self.__get_gen_plane(jl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def _get_cube(self):
        ''' Get Stokes profiles memmap
        '''

        # If not complete
        if not self.__complete:
            _error('The memmap can only be created if the file ' + \
                   'is complete')
            return None
        return self.__stk[:,:,:,::-1]

######################################################################
######################################################################
######################################################################

class hanlertio_class():
    ''' Class to manage files in the input or output of the HanleRT
        code
    '''

    def __init__(self,filename,silent=False):
        ''' Class initializer, process the filename
        '''

        # Arguments
        self.__verbosity = not silent

        # Check that it is a file
        if not os.path.isfile(filename):
            _error(f'{filename} is not a file',1,True)
            self.__valid = False
            return None

        # Save filename
        self.__filename = filename

        # Valid file?
        self.__valid = self.__get_class()

        # Not valid
        if not self.__valid:
            return None

    def valid(self):
        ''' Return if valid file
        '''
        return self.__valid

    def __verbose(self,msg,force=False):
        ''' Control verbosity
        '''
        if self.__verbosity or force: _verbose(msg)

    def __is_ascii(self):
        ''' Tries to identify if it is an ASCII file in
            one of the valid HanleRT formats
        '''

        # Open file
        f = open(self.__filename,'rb')

        # Read a bunch of bytes
        if 4096 < os.path.getsize(self.__filename):
            chunk = f.read(4096)
        else:
            chunk = f.read(os.path.getsize(self.__filename))

        # Close file
        f.close()

        # Try to interpret
        try:

            # First we check the first chunk to discard pure
            # binaries
            chunk = chunk.decode()

            # Initialize
            could_be_atmo = False
            could_be_atmo_b = False

            # At least partly ASCII, could still be a fits
            # and checked the header. We check line by line
            # now
            with open(self.__filename) as f:

                # Number of valid lines counter
                ival = 0

                # For each line
                for line in f:

                    # Strip comments
                    if '!' in line:
                        il = line.find('!')
                        cline = line[:il]
                    elif '*' in line:
                        il = line.find('*')
                        cline = line[:il]
                    else:
                        cline = line

                    # Strip spaces
                    cline = cline.strip()

                    # If empty, skip
                    if len(cline) < 1: continue

                    # There is content
                    ival += 1

                    # If first valid, check if 24 entries
                    if ival == 1:
                        if len(cline.split()) == 24:
                            could_be_atmo_b = True

                    # If more than one valid and could be
                    # second atmospheric ASCII format
                    if ival > 1 and could_be_atmo_b:
                        if len(cline.split()) != 24:
                            could_be_atmo_b = False

                    # If it is an atmospheric file, the
                    # second valid line should contain height
                    # or tau
                    if ival == 2:
                        if 'height' in cline.lower() or \
                           'tau' in cline.lower():
                            could_be_atmo = True
                        else:
                            could_be_atmo = False

                    # If it is an atmospheric file, the third
                    # valid line should contain a number
                    if ival == 3 and could_be_atmo:
                        try:
                            num = float(cline)
                        except:
                            could_be_atmo = False
                            raise

                    # If it is an atmospheric file, the fourth
                    # valid line should contain an integer
                    if ival == 4 and could_be_atmo:
                        try:
                            nz = int(cline)
                        except:
                            could_be_atmo = False
                            raise

                    # If it is an atmospheric file, the fifth
                    # valid line should contain at least
                    # five number
                    if ival >= 5 and ival < 5+nz and could_be_atmo:
                        try:
                            cols = cline.split()
                            if len(cols) >= 5:
                                for c in cols:
                                    try:
                                        num = float(c)
                                    except:
                                        could_be_atmo = False
                                        raise
                                        break
                            else:
                                could_be_atmo = False
                        except:
                            raise
                            could_be_atmo = False

                    # If valid atmospheric file, after the block
                    # we need to find another block or a correct
                    # label
                    if ival > 5 and could_be_atmo:
                        if ival == 5+nz:
                            try:
                                cols = cline.split()
                                nh = False
                                if len(cols) == 6:
                                    nh = True
                                    for c in cols:
                                        try:
                                            num = float(c)
                                        except:
                                            nh = False
                                            break
                                if not nh:
                                    if cline.lower() != 'ne' and \
                                       cline.lower() != 'rhoe' and \
                                       cline.lower() != 'pg' and \
                                       cline.lower() != 'pe' and \
                                       cline.lower() != 'rho':
                                        could_be_atmo = False
                            except:
                                could_be_atmo = False

                    # If it is an atmospheric file, the second
                    # block should have valid number if it
                    # exists
                    if ival > 5 and could_be_atmo:
                        if ival >= nz+5 and ival < 5+2*nz:
                            if nh:
                                try:
                                    cols = cline.split()
                                    if len(cols) == 6:
                                        for c in cols:
                                            try:
                                                num = float(c)
                                            except:
                                                could_be_atmo = False
                                                raise
                                                break
                                    else:
                                        could_be_atmo = False
                                except:
                                    raise
                                    could_be_atmo = False
                                        

            # If could be atmosphere
            if could_be_atmo:

                # If we read at least nz + 5 lines or 2*nz + 4 lines
                if ival == nz+5 or ival == 2*nz + 4:

                    # This is an atmosphere
                    return 0

            # If could be atmosphere second format
            if could_be_atmo_b:

                # This is an atmosphere in second format
                return 1

        except:
            raise
            return -1

    def __get_class(self):
        ''' Identify the type of file and link the suitable class
        '''

        # TODO TODO
        '''
  'CLEt': 'Optical depth from CLE synthesis', \
  'MRC': 'Maximum relative change from 1.5D synthesis', \
  'sp': 'Solution file with polarization', \
  'si': 'Solution file without polarization', \
        '''
        # TODO TODO

        # Possible labels
        labels = { \
          'invo': 'Inversion Result file', \
          'invi': 'Inversion input file', \
          '2Dbe': 'Emergent Stokes parameters from 1.5D synthesis', \
          '2Dbc': 'Contribution function from 1.5D synthesis', \
          '2Dbt': 'Height for optical depth unity from 1.5D ' + \
                  'synthesis', \
          '2Dct': 'Term to term collisional rates from 1.5D ' + \
                  'synthesis', \
          '2Dcl': 'Level to level collisional rates from 1.5D ' + \
                  'synthesis', \
          '2Dda': 'Damping parameters from 1.5D synthesis', \
          '2Dqe': 'Elastic rates from 1.5D synthesis', \
          '2Dba': 'Background continuum from 1.5D synthesis', \
          '2Dbp': 'Atomic populations from 1.5D synthesis', \
          '2Dbb': 'Departure coefficients from 1.5D synthesis', \
          '2Dat': 'multi-dimensional model atmosphere', \
          'CLEe': 'Emergent Stokes parameters from CLE synthesis', \
          'CLEt': 'Optical depth from CLE synthesis', \
          'MRC': 'Maximum relative change from 1.5D synthesis', \
          'sp': 'Solution file with polarization', \
          'si': 'Solution file without polarization', \
          'bj': 'Radiation field tensors from 1D synthesis', \
          'br': 'Density matrices from 1D synthesis', \
          'bp': 'Atomic populations from 1D synthesis', \
          'bb': 'Departure coefficients from 1D synthesis', \
          'bo': 'Stokes parameters in the quadrature in 1D ' + \
                'synthesis', \
          'ko': 'Frequency dependent radiation field tensors ' + \
                'in 1D synthesis', \
          'be': 'Emergent Stokes parameters from 1D synthesis', \
          'bc': 'Contribution function from 1D synthesis', \
          'bt': 'Height for optical depth unity from 1D synthesis', \
          'ct': 'Term to term collisional rates from 1D synthesis', \
          'cl': 'Level to level collisional rates from 1D ' + \
                'synthesis', \
          'da': 'Damping parameters from 1D synthesis', \
          'qe': 'Elastic rates from 1D synthesis', \
          'ba': 'Background continuum from 1D synthesis'}

        # Open the file and get two characters
        f = open(self.__filename,'rb')
        try:
          label = f.read(2).decode('utf-8')
        except:
          raise

        # Check if label within 2 character labels
        if label in labels:

          # Message
          self.__verbose(labels[label])

          # stokes_1D
          if label == 'be':

            # Close
            f.close()

            # Load Stokes 1D class
            self.__object = _stokes_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nl = self.__get_nl
              self.get_th = self.__get_th
              self.get_ph = self.__get_ph
              self.get_mu = self.__get_mu
              self.get_lambda = self.__get_lambda
              self.get_stokesi = self.__get_stokesi1d
              self.get_stokesq = self.__get_stokesq1d
              self.get_stokesu = self.__get_stokesu1d
              self.get_stokesv = self.__get_stokesv1d
              self.get_stokes = self.__get_stokes1d
              self.get_linear = self.__get_linear1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # stokesquad_1D
          elif label == 'bo':

            # Close
            f.close()

            # Load Stokes 1D class
            self.__object = _stokesquad_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nl = self.__get_nl
              self.get_nd = self.__get_nd
              self.get_th = self.__get_th
              self.get_ph = self.__get_ph
              self.get_mu = self.__get_mu
              self.get_lambda = self.__get_lambda
              self.get_stokesi = self.__get_stokesi1d
              self.get_stokesq = self.__get_stokesq1d
              self.get_stokesu = self.__get_stokesu1d
              self.get_stokesv = self.__get_stokesv1d
              self.get_stokes = self.__get_stokes1d
              self.get_linear = self.__get_linear1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # contribution_1D
          elif label == 'bc':

            # Close
            f.close()

            # Load contribution 1D class
            self.__object = _contribution_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nl = self.__get_nl
              self.get_nz = self.__get_nz
              self.get_th = self.__get_th
              self.get_ph = self.__get_ph
              self.get_mu = self.__get_mu
              self.get_lambda = self.__get_lambda
              self.get_height = self.__get_height
              self.get_ctri = self.__get_ctri1d
              self.get_ctrq = self.__get_ctrq1d
              self.get_ctru = self.__get_ctru1d
              self.get_ctrv = self.__get_ctrv1d
              self.get_ctr = self.__get_ctr1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # tau1_1D
          elif label == 'bt':

            # Close
            f.close()

            # Load tau1 1D class
            self.__object = _tau_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nl = self.__get_nl
              self.get_th = self.__get_th
              self.get_ph = self.__get_ph
              self.get_mu = self.__get_mu
              self.get_lambda = self.__get_lambda
              self.get_height = self.__get_height1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # jkq_1D
          elif label == 'bj':

            # Close
            f.close()

            # Load jkq 1D class
            self.__object = _jkq_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nz = self.__get_nz
              self.get_na = self.__get_na
              self.get_nt = self.__get_nt
              self.get_height = self.__get_height
              self.get_jkq = self.__get_jkq1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # rkq_1D
          elif label == 'br':

            # Close
            f.close()

            # Load rkq 1D class
            self.__object = _rkq_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nz = self.__get_nz
              self.get_na = self.__get_na
              self.get_height = self.__get_height
              self.get_rkq = self.__get_rkq1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # jkqnu_1D
          elif label == 'ko':

            # Close
            f.close()

            # Load jkq 1D class
            self.__object = _jkqnu_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nz = self.__get_nz
              self.get_nl = self.__get_nl
              self.get_lambda = self.__get_lambda
              self.get_jkq = self.__get_jkqnu1d

              # Valid class
              return True

            # Fail
            else:

              # Not valid class
              return False

          # Background
          elif label == 'ba':

            # Close
            f.close()

            # Load background 1D class
            self.__object = _back_1D(self.__filename)

            if self.__object is not None:

              # Methods
              self.get_filename = self.__get_filename
              self.get_nl = self.__get_nl
              self.get_nd = self.__get_nd
              self.get_nz = self.__get_nz
              self.get_dims = self.__get_dims
              self.get_vars = self.__get_vars
              self.get_vars_alias = self.__get_vars_alias
              self.get_vars_units = self.__get_vars_units
              self.get_lambda = self.__get_lambda
              self.get_data = self.__get_column_1dback

              # Valid class
              return True

          # Collisions, damping, and elastic rates
          elif label == 'ct' or label == 'cl' or \
               label == 'da' or label == 'qe':

            # Close
            f.close()

            # Load populations and departure 1.5D class
            self.__object = _cols_damp_1D(self.__filename)

            if self.__object is not None:

                # Methods
                self.get_filename = self.__get_filename
                self.get_type = self.__get_type
                self.get_nl = self.__get_nl
                self.get_nt = self.__get_nt
                self.get_nz = self.__get_nz
                self.get_dims = self.__get_dims
                self.get_data = self.__get_column_1dcdq

                # Valid class
                return True

          # Populations and departures
          elif label == 'bp' or label == 'bb':

            # Close
            f.close()

            # Load populations and departure 1.5D class
            self.__object = _pop_dep_1D(self.__filename)

            if self.__object is not None:

                # Methods
                self.get_filename = self.__get_filename
                self.get_type = self.__get_type
                self.get_nl = self.__get_nl
                self.get_nz = self.__get_nz
                self.get_dims = self.__get_dims
                self.get_data = self.__get_column_1dcapd

                # Valid class
                return True

          else:

            print(f' # {label} found in 2 char')

        # Not in two character labels
        else:

          # get one more character
          try:
            label += f.read(1).decode('utf-8')
          except:
            raise

          # Check if label within 3 character labels
          if label in labels:

            # Message
            self.__verbose(labels[label])

            print(f' # {label} found in 3 char')

          # Not in three character labels
          else:

            # get one more character
            try:
              label += f.read(1).decode('utf-8')
            except:
              raise

            # Check if label within 4 character labels
            if label in labels:

              # Message
              self.__verbose(labels[label])

              # 1.5D Stokes
              if label == '2Dbe':

                # Close
                f.close()

                # Load Stokes 1.5D class
                self.__object = _stokes_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_polarized = self.__get_polarized
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nl = self.__get_nl
                  self.get_nxy = self.__get_nxy
                  self.get_dims = self.__get_dims
                  self.get_th = self.__get_th
                  self.get_ph = self.__get_ph
                  self.get_mu = self.__get_mu
                  self.get_lambda = self.__get_lambda
                  self.get_stokesi_column = self.__get_stokesi15d_c
                  self.get_stokesq_column = self.__get_stokesq15d_c
                  self.get_stokesu_column = self.__get_stokesu15d_c
                  self.get_stokesv_column = self.__get_stokesv15d_c
                  self.get_stokes_column = self.__get_stokes15d_c
                  self.get_linear_column = self.__get_linear15d_c
                  self.get_stokesi_plane = self.__get_stokesi15d_p
                  self.get_stokesq_plane = self.__get_stokesq15d_p
                  self.get_stokesu_plane = self.__get_stokesu15d_p
                  self.get_stokesv_plane = self.__get_stokesv15d_p
                  self.get_stokes_plane = self.__get_stokes15d_p
                  self.get_linear_plane = self.__get_linear15d_p
                  self.get_scale = self.__get_scale
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # 1.5D Contribution function
              elif label == '2Dbc':

                # Close
                f.close()

                # Load contribution 1.5D class
                self.__object = _contribution_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nz = self.__get_nz
                  self.get_nl = self.__get_nl
                  self.get_nxy = self.__get_nxy
                  self.get_nxyz = self.__get_nxyz
                  self.get_dims = self.__get_dims
                  self.get_th = self.__get_th
                  self.get_ph = self.__get_ph
                  self.get_mu = self.__get_mu
                  self.get_lambda = self.__get_lambda
                  self.get_ctri_column = self.__get_ctri15d_c
                  self.get_ctrq_column = self.__get_ctrq15d_c
                  self.get_ctru_column = self.__get_ctru15d_c
                  self.get_ctrv_column = self.__get_ctrv15d_c
                  self.get_ctr_column = self.__get_ctr15d_c
                  self.get_scale = self.__get_scale
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # 1.5D height tau equal 1
              elif label == '2Dbt':

                # Close
                f.close()

                # Load Tau 1.5D class
                self.__object = _tau_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nxy = self.__get_nxy
                  self.get_nl = self.__get_nl
                  self.get_dims = self.__get_dims
                  self.get_th = self.__get_th
                  self.get_ph = self.__get_ph
                  self.get_mu = self.__get_mu
                  self.get_lambda = self.__get_lambda
                  self.get_column = self.__get_height15d_c
                  self.get_plane = self.__get_height15d_p
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # Model atmosphere
              elif label == '2Dat':

                # Close
                f.close()

                # Load 3D atmospheric model class
                self.__object = _atmo_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_precision = self.__get_precision
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nz = self.__get_nz
                  self.get_nxyz = self.__get_nxyz
                  self.get_dims = self.__get_dims
                  self.get_vars = self.__get_vars
                  self.get_vars_units = self.__get_vars_units
                  self.get_vars_alias = self.__get_vars_alias
                  self.get_column = self.__get_column_2dat
                  self.get_plane = self.__get_plane_2dat
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # inversion_in
              elif label == 'invi':

                # Close
                f.close()

                # Load inversion output class
                self.__object = _inversion_in(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_info_verb = self.__get_info_verb
                  self.get_info = self.__get_info
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nxy = self.__get_nxy
                  self.get_nl = self.__get_nl
                  self.get_dims = self.__get_dims
                  self.get_lambda = self.__get_lambda
                  self.get_los = self.__get_los
                  self.get_los_column = self.__get_los_c
                  self.get_los_plane = self.__get_los_p
                  self.get_stokesi_column = self.__get_stokesi15d_c
                  self.get_stokesq_column = self.__get_stokesq15d_c
                  self.get_stokesu_column = self.__get_stokesu15d_c
                  self.get_stokesv_column = self.__get_stokesv15d_c
                  self.get_linear_column = self.__get_linear15d_c
                  self.get_stokes_column = self.__get_stokes15d_c
                  self.get_stokesi_plane = self.__get_stokesi15d_p
                  self.get_stokesq_plane = self.__get_stokesq15d_p
                  self.get_stokesu_plane = self.__get_stokesu15d_p
                  self.get_stokesv_plane = self.__get_stokesv15d_p
                  self.get_linear_plane = self.__get_linear15d_p
                  self.get_stokes_plane = self.__get_stokes15d_p
                  self.get_sigma = self.__get_sigma
                  self.get_sigmai_column = self.__get_sigmai15d_c
                  self.get_sigmaq_column = self.__get_sigmaq15d_c
                  self.get_sigmau_column = self.__get_sigmau15d_c
                  self.get_sigmav_column = self.__get_sigmav15d_c
                  self.get_sigma_column = self.__get_sigma15d_c
                  self.get_sigmai_plane = self.__get_sigmai15d_p
                  self.get_sigmaq_plane = self.__get_sigmaq15d_p
                  self.get_sigmau_plane = self.__get_sigmau15d_p
                  self.get_sigmav_plane = self.__get_sigmav15d_p
                  self.get_sigma_plane = self.__get_sigma15d_p
                  self.get_diffi_column = self.__get_diffi15d_c
                  self.get_diffq_column = self.__get_diffq15d_c
                  self.get_diffu_column = self.__get_diffu15d_c
                  self.get_diffv_column = self.__get_diffv15d_c
                  self.get_diff_column = self.__get_diff15d_c
                  self.get_diffi_plane = self.__get_diffi15d_p
                  self.get_diffq_plane = self.__get_diffq15d_p
                  self.get_diffu_plane = self.__get_diffu15d_p
                  self.get_diffv_plane = self.__get_diffv15d_p
                  self.get_diff_plane = self.__get_diff15d_p
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # inversion_out
              elif label == 'invo':

                # Close
                f.close()

                # Load inversion output class
                self.__object = _inversion_out(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_polarized = self.__get_polarized
                  self.get_vtype = self.__get_vtype
                  self.get_btype = self.__get_btype
                  self.get_jkqin = self.__get_jkqin
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nz = self.__get_nz
                  self.get_nl = self.__get_nl
                  self.get_nxy = self.__get_nxy
                  self.get_nxyz = self.__get_nxyz
                  self.get_dims = self.__get_dims
                  self.get_nvar_atmo = self.__get_nvar_atmo
                  self.get_nvar = self.__get_nvar
                  self.get_vars = self.__get_vars
                  self.get_vars_units = self.__get_vars_units
                  self.get_vars_fix = self.__get_vars_fix
                  self.get_vars_fix_units = self.__get_vars_fix_units
                  self.get_vars_atmo = self.__get_vars_atmo
                  self.get_vars_atmo_units= self.__get_vars_atmo_units
                  self.get_lambda = self.__get_lambda
                  self.get_column = self.__get_column_res
                  self.get_column_atmo = self.__get_column_atmo
                  self.get_column_rf = self.__get_column_rf
                  self.get_plane_chi = self.__get_plane_chi
                  self.get_plane_atmo = self.__get_plane_atmo
                  self.get_plane_stk = self.__get_plane_stk
                  self.get_node = self.__get_node
                  self.get_cube = self.__get_cube
                  self.get_cube_atmo = self.__get_cube_atmo

                  # Valid class
                  return True

              # collisional rates, damping, and elastic rates
              elif label == '2Dct' or label == '2Dcl' or \
                   label == '2Dda' or label == '2Deq':

                # Close
                f.close()

                # Load collisions and damping 1.5D class
                self.__object = _cols_damp_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  if label == '2Dct' or label == '2Dcl':
                      self.get_type = self.__get_type
                  self.get_nentry = self.__get_nentry
                  self.get_na = self.__get_na
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nz = self.__get_nz
                  self.get_nxy = self.__get_nxy
                  self.get_nxyz = self.__get_nxyz
                  self.get_dims = self.__get_dims
                  self.get_column = self.__get_column_15dcapd
                  self.get_plane = self.__get_plane_15dcapd
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # Continuum quantities
              elif label == '2Dba':

                # Close
                f.close()

                # Load continuum 1.5D class
                self.__object = _back_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nz = self.__get_nz
                  self.get_nl = self.__get_nl
                  self.get_nxy = self.__get_nxy
                  self.get_nxyz = self.__get_nxyz
                  self.get_dims = self.__get_dims
                  self.get_vars = self.__get_vars
                  self.get_vars_units = self.__get_vars_units
                  self.get_vars_alias= self.__get_vars_alias
                  self.get_lambda = self.__get_lambda
                  self.get_column = self.__get_column_res
                  self.get_plane = self.__get_plane_back
                  self.get_scale = self.__get_scale
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # Populations and departures
              elif label == '2Dbp' or label == '2Dbb':

                # Close
                f.close()

                # Load populations and departure 1.5D class
                self.__object = _pop_dep_15D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_type = self.__get_type
                  self.get_nentry = self.__get_nentry
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nz = self.__get_nz
                  self.get_nxy = self.__get_nxy
                  self.get_nxyz = self.__get_nxyz
                  self.get_dims = self.__get_dims
                  self.get_column = self.__get_column_15dcapd
                  self.get_plane = self.__get_plane_15dcapd
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # CLE Stokes
              elif label == 'CLEe':

                # Close
                f.close()

                # Load populations and departure 1.5D class
                self.__object = _stokes_CLE(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_type = self.__get_type
                  self.get_nx = self.__get_nx
                  self.get_ny = self.__get_ny
                  self.get_nxy = self.__get_nxy
                  self.get_nl = self.__get_nl
                  self.get_dims = self.__get_dims
                  self.get_lambda = self.__get_lambda
                  self.get_geometry = self.__get_geometry
                  self.get_stokesi_column = self.__get_stokesi15d_c
                  self.get_stokesq_column = self.__get_stokesq15d_c
                  self.get_stokesu_column = self.__get_stokesu15d_c
                  self.get_stokesv_column = self.__get_stokesv15d_c
                  self.get_stokes_column = self.__get_stokes15d_c
                  self.get_linear_column = self.__get_linear15d_c
                  self.get_stokesi_plane = self.__get_stokesi15d_p
                  self.get_stokesq_plane = self.__get_stokesq15d_p
                  self.get_stokesu_plane = self.__get_stokesu15d_p
                  self.get_stokesv_plane = self.__get_stokesv15d_p
                  self.get_stokes_plane = self.__get_stokes15d_p
                  self.get_linear_plane = self.__get_linear15d_p
                  self.get_scale = self.__get_scale
                  self.get_cube = self.__get_cube

                  # Valid class
                  return True

              # Fail
              else:

                # Not valid class
                return False

                print(f' # {label} found in 4 char')

            # Not any of the labels
            else:

              # Close
              f.close()

              # Try to identify an ASCII format
              label = self.__is_ascii()

              #
              # Normal 1D model atmosphere
              if label == 0:

                # Load 3D atmospheric model class
                self.__object = _atmo_1D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_name = self.__get_name
                  self.get_comment = self.__get_comment
                  self.get_nz = self.__get_nz
                  self.get_wavelength = self.__get_wavelength
                  self.get_vars = self.__get_vars
                  self.get_vars_units = self.__get_vars_units
                  self.get_vars_alias = self.__get_vars_alias
                  self.get_column = self.__get_column_1dat

                  # Valid class
                  return True

              #
              # Full 1D model atmosphere
              elif label == 1:

                # Load 3D atmospheric model class
                self.__object = _atmo_b_1D(self.__filename)

                if self.__object is not None:

                  # Methods
                  self.get_filename = self.__get_filename
                  self.get_nz = self.__get_nz
                  self.get_vars = self.__get_vars
                  self.get_vars_units = self.__get_vars_units
                  self.get_vars_alias = self.__get_vars_alias
                  self.get_column = self.__get_column_1dat

                  # Valid class
                  return True


        # Return empty
        return False

    def help(self,method=None,input=None):
        ''' Get help for the methods and input parameters for the
            stablished class
        '''
        if not self.__valid:
          self.__verbose('No valid file has been provided',True)
        else:

          # Get data about methods
          methods = self.__object._get_help()

          # If nothing specified, just dump methods
          if method is None or method not in methods:

            # Dump method names
            self.__verbose('Available methods:',True)
            for m in methods:
                self.__verbose(f'  - {m}',True)
            self.__verbose('You can call ' + \
                           'help(method=name,input=name) ' + \
                           'to get information about a ' + \
                           'method or input',True)

          # If methods specified, provide particular info
          else:

            # Decide branch
            noin = input is None or methods[method][0] is None
            if not noin:
                noin = input not in methods[method][0]

            # If no input specified, give info
            if noin:

              # Print help
              self.__verbose(f'{method}: {methods[method][1]}',True)

              # Has inputs?
              if methods[method][0] is not None:
                  self.__verbose(f'Inputs:',True)
                  for i in methods[method][0]:
                      self.__verbose(f'  - {i}',True)

            # Specified an input
            else:

              # Print help
              self.__verbose(f'In method {method}, ' + \
                             f'{input}: ' + \
                             f'{methods[method][0][input]}',True)

######################################################################

    # Parsers

    # All
    def __get_filename(self):
        return self.__object._get_filename()

    # 1D atmo
    def __get_name(self):
        return self.__object._get_name()

    # 1D atmo
    def __get_comment(self):
        return self.__object._get_comment()

    # 1D atmo
    def __get_wavelength(self):
        return self.__object._get_wavelegnth()

    # 15D atmo
    def __get_precision(self):
        return self.__object._get_precision()

    # 1D cols, 1D popdep, 15D cols, 15D popdep
    def __get_type(self):
        return self.__object._get_type()

    # stokes 15D, inversion out
    def __get_polarized(self):
        return self.__object._get_polarized()

    # inversion out
    def __get_vtype(self):
        return self.__object._get_vtype()

    # inversion out
    def __get_btype(self):
        return self.__object._get_btype()

    # inversion out
    def __get_jkqin(self):
        return self.__object._get_jkqin()

    # inversion in
    def __get_info_verb(self):
        return self.__object._get_info_verb()

    # inversion in
    def __get_info(self):
        return self.__object._get_info()

    # inversion out
    def __get_vars_fix(self):
        return self.__object._get_vars_fix()

    # inversion out
    def __get_vars_fix_units(self):
        return self.__object._get_vars_fix_units()

    # 1D atmos, 3D atmos, inversion out, 15D back, 1D back
    def __get_vars(self):
        return self.__object._get_vars()

    # 1D atmos, 3D atmos, inversion out, 15D back, 1D back
    def __get_vars_units(self):
        return self.__object._get_vars_units()

    # 1D atmos, 3D atmos, 15D back
    def __get_vars_alias(self):
        return self.__object._get_vars_alias()

    # inversion out
    def __get_vars_atmo(self):
        return self.__object._get_vars_atmo()

    # inversion out
    def __get_vars_atmo_units(self):
        return self.__object._get_vars_atmo_units()

    # stokes 1D, stokesquad 1D, contribution 1D, tau 1D, 1D back,
    # 1D popdep, jkqnu 1D, stokes 15D, contribution 15D, tau15D,
    # inversion in, inversion out, 15D back
    def __get_nl(self):
        return self.__object._get_nl()

    # 1D back, stokesquad 1D
    def __get_nd(self):
        return self.__object._get_nd()

    # stokes 15D, contribution 15D, tau15D, inversion in,
    # inversion out, 15D cols, 15D back, 15D popdep
    def __get_nx(self):
        return self.__object._get_nx()

    # stokes 15D, contribution 15D, tau15D, inversion in,
    # inversion out, 15D cols, 15D back, 15D popdep
    def __get_ny(self):
        return self.__object._get_ny()

    # stokes 15D, contribution 15D, tau15D, inversion in,
    # inversion out, 15D cols, 15D back, 15D popdep
    def __get_nxy(self):
        return self.__object._get_nxy()

    # contribution 15D, inversion out, 15D cols, 15D back, 15D popdep
    def __get_nxyz(self):
        return self.__object._get_nxyz()

    # stokes 15D, contribution 15D, tau15D, 3D atmos, inversion in,
    # inversion out, 15D cols, 15D back, 15D popdep, 1D cols, 1D back,
    # 1D popdep
    def __get_dims(self):
        return self.__object._get_dims()

    # Stokes 15D, contribution 15D, 15D back, CLE
    def __get_scale(self):
        return self.__object._get_scale()

    # contribution 1D, jkq 1D, rkq 1D, jkqnu 1D, 1D popdep,
    # contribution 15D, 1D atmos, 3D atmos, inversion out, 15D cols,
    # 15D back, 15D popdep
    def __get_nz(self):
        return self.__object._get_nz()

    # 15D cols, 15D popdep
    def __get_nentry(self):
        return self.__object._get_nentry()

    # jkq 1D, rkq 1D, 15D cols
    def __get_na(self):
        return self.__object._get_na()

    # jkq 1D, 1D damp
    def __get_nt(self,atom=None):
        return self.__object._get_nt(atom)

    # inversion out
    def __get_nvar_atmo(self):
        return self.__object._get_nvar_atmo()

    # inversion out
    def __get_nvar(self):
        return self.__object._get_nvar()

    # inversion in
    def __get_los(self):
        return self.__object._get_los()

    # inversion in
    def __get_los_c(self,ix=None,iy=None):
        return self.__object._get_los_column(ix,iy)

    # inversion in
    def __get_los_p(self,ix=None,iy=None):
        return self.__object._get_los_plane()

    # stokes 1D, stokesquad 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D
    def __get_th(self):
        return self.__object._get_th()

    # stokes 1D, stokesquad 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D
    def __get_ph(self):
        return self.__object._get_ph()

    # stokes 1D, stokesquad 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D
    def __get_mu(self):
        return self.__object._get_mu()

    # stokes 1D, stokesquad 1D, contribution 1D, tau 1D, jkqnu 1D,
    # stokes 15D, contribution 15D, tau15D, inversion in, inversion
    # out 15D back, 1D back
    def __get_lambda(self,minl=None,maxl=None):
        return self.__object._get_lambda(minl,maxl)

    # CLEe
    def __get_geometry(self):
        return self.__object._get_geometry()

    # contribution 1D, jkq 1D, rkq 1D
    def __get_height(self,minh=None,maxh=None):
        return self.__object._get_height(minh,maxh)

    # tau 1D
    def __get_height1d(self,minl=None,maxl=None):
        return self.__object._get_height(minl,maxl)

    # stokes 1D, stokesquad 1D
    def __get_stokesi1d(self,minl=None,maxl=None):
        return self.__object._get_stokesi(minl,maxl)

    # stokes 1D, stokesquad 1D
    def __get_stokesq1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesq(minl,maxl,fractional)

    # stokes 1D, stokesquad 1D
    def __get_stokesu1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesu(minl,maxl,fractional)

    # stokes 1D, stokesquad 1D
    def __get_stokesv1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesv(minl,maxl,fractional)

    # stokes 1D, stokesquad 1D
    def __get_stokes1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokes(minl,maxl,fractional)

    # stokes 1D, stokesquad 1D
    def __get_linear1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_linear(minl,maxl,fractional)

    # contribution 1D
    def __get_ctri1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctri(minl,maxl,minh,maxh)

    # contribution 1D
    def __get_ctrq1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctrq(minl,maxl,minh,maxh)

    # contribution 1D
    def __get_ctru1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctru(minl,maxl,minh,maxh)

    # contribution 1D
    def __get_ctrv1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctrv(minl,maxl,minh,maxh)

    # contribution 1D
    def __get_ctr1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctr(minl,maxl,minh,maxh)

    # jkq 1D
    def __get_jkq1d(self,atom=None,transition=None,k=None,q=None, \
                    minh=None,maxh=None,sti=False):
        return self.__object._get_jkq(atom,transition,k,q, \
                                      minh,maxh,sti)

    # jkqnu 1D
    def __get_jkqnu1d(self,minl=None,maxl=None):
        return self.__object._get_jkq(minl,maxl)

    # rkq 1D
    def __get_rkq1d(self,atom=None,minh=None,maxh=None):
        return self.__object._get_rkq(atom,minh,maxh)

    # contribution 1D
    def __get_ctr1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctr(minl,maxl,minh,maxh)

    # 1D cols, damping, qel
    def __get_column_1dcdq(self,ia=None,i1=None,i2=None,it=None):
        return self.__object._get_data(ia,i1,i2,it)

    # 1D back
    def __get_column_1dback(self,minl=None,maxl=None,var=None):
        return self.__object._get_data(minl,maxl,var)

    # 1D popdep
    def __get_column_1dcapd(self,ie=None):
        return self.__object._get_data(ie)

    # stokes 15D, inversion in
    def __get_stokesi15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_stokesi_column(ix,iy,minl,maxl)

    # stokes 15D, inversion in
    def __get_stokesq15d_c(self,ix=None,iy=None, \
                           minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesq_column(ix,iy,minl,maxl, \
                                                 fractional)

    # stokes 15D, inversion in
    def __get_stokesu15d_c(self,ix=None,iy=None, \
                           minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesu_column(ix,iy,minl,maxl, \
                                                 fractional)

    # stokes 15D, inversion in
    def __get_stokesv15d_c(self,ix=None,iy=None, \
                           minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesv_column(ix,iy,minl,maxl, \
                                                 fractional)

    # stokes 15D, inversion in
    def __get_stokes15d_c(self,ix=None,iy=None, \
                          minl=None,maxl=None,fractional=False):
        return self.__object._get_stokes_column(ix,iy,minl,maxl, \
                                                fractional)

    # stokes 15D, inversion in
    def __get_linear15d_c(self,ix=None,iy=None, \
                          minl=None,maxl=None,fractional=False):
        return self.__object._get_linear_column(ix,iy,minl,maxl, \
                                                fractional)

    # stokes 15D, inversion in
    def __get_stokesi15d_p(self,il=None):
        return self.__object._get_stokesi_plane(il)

    # stokes 15D, inversion in
    def __get_stokesq15d_p(self,il=None,fractional=False):
        return self.__object._get_stokesq_plane(il,fractional)

    # stokes 15D, inversion in
    def __get_stokesu15d_p(self,il=None,fractional=False):
        return self.__object._get_stokesu_plane(il,fractional)

    # stokes 15D, inversion in
    def __get_stokesv15d_p(self,il=None,fractional=False):
        return self.__object._get_stokesv_plane(il,fractional)

    # stokes 15D, inversion in
    def __get_stokes15d_p(self,il=None,fractional=False):
        return self.__object._get_stokes_plane(il,fractional)

    # stokes 15D, inversion in
    def __get_linear15d_p(self,il=None,fractional=False):
        return self.__object._get_linear_plane(il,fractional)

    # contribution 15D
    def __get_ctri15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_ctri_column(ix,iy,minl,maxl)

    # contribution 15D
    def __get_ctrq15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_ctrq_column(ix,iy,minl,maxl)

    # contribution 15D
    def __get_ctru15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_ctru_column(ix,iy,minl,maxl)

    # contribution 15D
    def __get_ctrv15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_ctrv_column(ix,iy,minl,maxl)

    # contribution 15D
    def __get_ctr15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_ctr_column(ix,iy,minl,maxl)

    # Stokes 15D, contribution 15D, tau15D, 15D atmo, inversion in,
    # inversion out, 15D cols, 15D back, 15D popdep, CLE
    def __get_cube(self):
        return self.__object._get_cube()

    # inversion out
    def __get_cube_atmo(self):
        return self.__object._get_cube_atmo()

    # tau 15D
    def __get_height15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_column(ix,iy,minl,maxl)

    # tau 15D
    def __get_height15d_p(self,il=None):
        return self.__object._get_height_plane(il)

    # 1D atmo
    def __get_column_1dat(self,minz=None,maxz=None,var=None):
        return self.__object._get_column(minz,maxz,var)

    # 3D atmo
    def __get_column_2dat(self,ix=None,iy=None,minh=None,maxh=None, \
                          mint=None,maxt=None,var=None):
        return self.__object._get_column(ix,iy,minh,maxh, \
                                         mint,maxt,var)

    # 3D atmo
    def __get_plane_2dat(self,iz,var=None):
        return self.__object._get_plane(iz,var)

    # inversion out
    def __get_column_atmo(self,ix,iy,minh=None,maxh=None,var=None):
        return self.__object._get_column_atmo(ix,iy,minh,maxh,var)

    # inversion out, 15D back
    def __get_column_res(self,ix,iy,minl=None,maxl=None,var=None):
        return self.__object._get_column(ix,iy,minl,maxl,var)

    # inversion out
    def __get_column_rf(self,ix,iy,minl=None,maxl=None,var=None):
        return self.__object._get_column_rf(ix,iy,minl,maxl,var)

    # inversion out
    def __get_plane_chi(self):
        return self.__object._get_plane_chi()

    # inversion out
    def __get_plane_stk(self,il,var=None):
        return self.__object._get_plane_stk(il,var)

    # inversion out
    def __get_plane_atmo(self,iz,var=None):
        return self.__object._get_plane_atmo(iz,var)

    # inversion out
    def __get_node(self,var):
        return self.__object._get_node(var)

    # 15D back
    def __get_plane_back(self,iz,minl=None,maxl=None,var=None):
        return self.__object._get_plane(iz,minl,maxl,var)

    # 15D cols, 15D popdeb
    def __get_column_15dcapd(self,ix=None,iy=None,ie=None):
        return self.__object._get_column(ix,iy,ie)

    # 15D cols, 15D popbed
    def __get_plane_15dcapd(self,iz,ie=None):
        return self.__object._get_plane(iz,ie)

    # inversion in
    def __get_sigma(self):
        return self.__object._get_sigma()

    # inversion in
    def __get_sigmai15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_sigmai_column(ix,iy,minl,maxl)

    # inversion in
    def __get_sigmaq15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_sigmaq_column(ix,iy,minl,maxl)

    # inversion in
    def __get_sigmau15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_sigmau_column(ix,iy,minl,maxl)

    # inversion in
    def __get_sigmav15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_sigmav_column(ix,iy,minl,maxl)

    # inversion in
    def __get_sigma15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_sigma_column(ix,iy,minl,maxl)

    # inversion in
    def __get_sigmai15d_p(self,il=None):
        return self.__object._get_sigmai_plane(il)

    # inversion in
    def __get_sigmaq15d_p(self,il=None):
        return self.__object._get_sigmaq_plane(il)

    # inversion in
    def __get_sigmau15d_p(self,il=None):
        return self.__object._get_sigmau_plane(il)

    # inversion in
    def __get_sigmav15d_p(self,i=None):
        return self.__object._get_sigmav_plane(il)

    # inversion in
    def __get_sigma15d_p(self,il=None):
        return self.__object._get_sigma_plane(il)

    # inversion in
    def __get_diffi15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_diffi_column(ix,iy,minl,maxl)

    # inversion in
    def __get_diffq15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_diffq_column(ix,iy,minl,maxl)

    # inversion in
    def __get_diffu15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_diffu_column(ix,iy,minl,maxl)

    # inversion in
    def __get_diffv15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_diffv_column(ix,iy,minl,maxl)

    # inversion in
    def __get_diff15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_diff_column(ix,iy,minl,maxl)

    # inversion in
    def __get_diffi15d_p(self,il=None):
        return self.__object._get_diffi_plane(il)

    # inversion in
    def __get_diffq15d_p(self,il=None):
        return self.__object._get_diffq_plane(il)

    # inversion in
    def __get_diffu15d_p(self,il=None):
        return self.__object._get_diffu_plane(il)

    # inversion in
    def __get_diffv15d_p(self,i=None):
        return self.__object._get_diffv_plane(il)

    # inversion in
    def __get_diff15d_p(self,il=None):
        return self.__object._get_diff_plane(il)

######################################################################
######################################################################
######################################################################
