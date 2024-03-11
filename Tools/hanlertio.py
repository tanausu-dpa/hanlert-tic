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

################################################################################
################################################################################
################################################################################

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
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__th = struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)
            self.__hsize = 22
            f.close()
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
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
        except:
            raise

    def __get_gen_stokes(self,minl=None,maxl=None,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters
        '''

        # Initialize
        out = [None,None,None,None]
        bsiz = self.__nl*8

        try:

            # Open and seek data
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)

            # Get wavelength
            if fractional or minl is not None or maxl is not None:
                omg = np.array(struct.unpack('d'*self.__nl, \
                                             f.read(bsiz)))
                lam = 1e2/omg[::-1]
            else:
                f.seek(bsiz,1)

            # Intensity
            if 0 in indx or fractional:

                # Get intensity
                stkI = np.array(struct.unpack('d'*self.__nl, \
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

                # If in output
                if j in indx:

                    # Read Stokes
                    out[j] = np.array(struct.unpack('d'*self.__nl, \
                                                    f.read(bsiz)))[::-1]
                # Not in output
                else:

                    # Skip
                    f.seek(bsiz,1)

                # If fractional and output
                if fractional and j in indx:
                    out[j] /= stkI

            # Close unit
            f.close()

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

        except struct.error:
            raise
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

################################################################################
################################################################################
################################################################################

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
           'Get vertical axis, either heights in [km] or in optical depth'], \
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
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
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
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            f.seek(self.__nl*8,1)
            z = np.array(struct.unpack('d'*self.__nz, \
                                         f.read(8*self.__nz)))
            if iminh is not None:
                i = np.argmin(np.absolute(z - iminh))
                z = z[i:]
            if imaxh is not None:
                i = np.argmin(np.absolute(z - imaxh))
                z = z[:i+1]
            f.close()
            return z
        except struct.error:
            raise
        except:
            raise

    def __get_gen_ctr(self,minl=None,maxl=None,minh=None,maxh=None,indx=[0]):
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
        bsiz = self.__nl*8

        # Get data
        try:

            # Open and seek data
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)

            # Cutting lambda
            if minl is not None or maxl is not None:
                omg = np.array(struct.unpack('d'*self.__nl, \
                                             f.read(self.__nl*8)))
                lam = 1e2/omg[::-1]
            # No cut
            else:
                # Skip
                f.seek(self.__nl*8,1)

            # Cutting height
            if minh is not None or maxh is not None:
                z = np.array(struct.unpack('d'*self.__nz, \
                                             f.read(8*self.__nz)))
            # Not cutting height
            else:
                f.seek(self.__nz*8,1)

            # For each Stokes
            for j in range(4):

                # Read
                if j in indx:

                    out[j] = (np.array( \
                          struct.unpack('d'*self.__nl*self.__nz, \
                                        f.read(8*self.__nl*self.__nz))). \
                             reshape((self.__nz,self.__nl)))[:,::-1]

                # No read
                else:

                    # Skip
                    f.seek(self.__nl*self.__nz*8,1)

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

            # Close
            f.close()

            # Units
            for j in indx:
                out[j] *= self.__unit_trans

            # Return
            return out

        except struct.error:
            raise
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

################################################################################
################################################################################
################################################################################

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
        ''' Reads hanlert 1D height where optical depth is one file head
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nl = struct.unpack('i',f.read(4))[0]
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__th = struct.unpack('d',f.read(8))[0]
            self.__ph = struct.unpack('d',f.read(8))[0]
            self.__mu = np.cos(np.pi*self.__th/180e0)
            self.__hsize = 26
            f.close()
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
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
        except:
            raise

    def _get_height(self,minl=None,maxl=None):
        ''' Get height for optical depth equal to 1
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            tau1 = np.array(struct.unpack('d'*self.__nl, \
                                          f.read(8*self.__nl)))
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                tau1 = tau1[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                tau1 = tau1[:i+1]
            f.close()
            return tau1
        except struct.error:
            raise
        except:
            raise

################################################################################
################################################################################
################################################################################

class _jkq_1D():
    ''' Class to manage the frequency integrated radiation field tensors
        from a 1D synthesis
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
           'Get vertical axis, either heights in [km] or in optical ' + \
           'depth'], \
          'get_jkq': \
          [{'atom': 'Request this particular atom index', \
            'transition': 'Request this particular transition index. ' + \
                          'Requires "atom"', \
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
            # Z
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
            f.close()
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
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            z = np.array(struct.unpack('d'*self.__nz, \
                                         f.read(8*self.__nz)))
            if iminh is not None:
                i = np.argmin(np.absolute(z - iminh))
                z = z[i:]
            if imaxh is not None:
                i = np.argmin(np.absolute(z - imaxh))
                z = z[:i+1]
            f.close()
            return z
        except struct.error:
            raise
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
        try:

            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            z = np.array(struct.unpack('d'*self.__nz, \
                                         f.read(8*self.__nz)))

            # If asking for stimulated, skip normal
            if sti:
                f.seek(np.sum(self.__jump+4),1)

            # No atom specified
            if atom is None:

                # Big output
                jkq = np.zeros((np.sum(self.__nt),3,5,2,self.__nz))

                # Running index
                ind = -1

                # For each atom
                for ia in range(self.__na):

                    # Skip transition
                    if not sti:
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
                    if sti:
                        f.seek(np.sum(self.__jump[:atom]),1)
                    else:
                        f.seek(np.sum(self.__jump[:atom]+4),1)

                # Skip transition
                if not sti:
                    f.seek(4,1)

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
                            jq = q + 2 - k

                            # Trim
                            jkq = jkq[jq,:,:]
            f.close()
            return jkq*self.__unit_trans
        except struct.error:
            raise
        except:
            raise

################################################################################
################################################################################
################################################################################

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
           'Get vertical axis, either heights in [km] or in optical ' + \
           'depth'], \
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
            # Get actual header
            f = open(self.__filename,'rb')
            f.seek(2,0)
            self.__nz = struct.unpack('i',f.read(4))[0]
            self.__na = struct.unpack('i',f.read(4))[0]
            # Z
            z = np.array(struct.unpack('d'*self.__nz, \
                                       f.read(8*self.__nz)))
            # Check if tau
            self.__ltau = np.min(z) > 0 and np.max(z) < 1e2
            self.__zreverse = z[-2] > z[-1]
            self.__hsize = 10
            f.close()
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
            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            z = np.array(struct.unpack('d'*self.__nz, \
                                       f.read(8*self.__nz)))
            if not self.__ltau: z *= 1e-5
            if iminh is not None:
                i = np.argmin(np.absolute(z - iminh))
                z = z[i:]
            if imaxh is not None:
                i = np.argmin(np.absolute(z - imaxh))
                z = z[:i+1]
            f.close()
            return z
        except struct.error:
            raise
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

        try:

            f = open(self.__filename,'rb')
            f.seek(self.__hsize,0)
            z = np.array(struct.unpack('d'*self.__nz, \
                                       f.read(8*self.__nz)))
            if not self.__ltau: z *= 1e-5

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
                    lrkq['n'] = np.array(struct.unpack('d'*self.__nz, \
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
                            if J1 not in lrkq[it]['J']: lrkq[it]['J'].append(J1)
                        J2 = struct.unpack('i',f.read(4))[0]
                        if keep:
                            if J2 not in lrkq[it]['J']: lrkq[it]['J'].append(J2)

                        # Get label
                        if keep:
                            if J1 % 2 == 0:
                                ijs1 = f'{int(round(J1*0.5))}'
                            else:
                                ijs1 = f'{J1}/2'
                            if ijs1 not in lrkq[it]: lrkq[it][ijs1] = {}
                            if mt:
                                if J2 % 2 == 0:
                                    ijs2 = f'{int(round(J2*0.5))}'
                                else:
                                    ijs2 = f'{J2}/2'
                                lrkq[it][ijs1][ijs2] = {}

                        # K limits
                        kmin = int(np.absolute(int(round(J1 - J2))//2))
                        kmax = int(np.absolute(int(round(J1 + J2))//2))

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
                                    point[K][Q] = np.zeros((self.__nz), \
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

################################################################################
################################################################################
################################################################################

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
          [None,' Get number of positions in height axis and levels'], \
         'get_data': \
          [{'ie': \
            'List of levels to include in the output'}, \
            'Extract the populations/departures']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert collisional/damping file head
        '''
        try:

            # Get actual header
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

            # Size of column
            self.__column = 8*self.__nz*self.__nn

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
        ''' Get number of positions in x, y, and height axes and entries
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

        # Size column
        size = self.__column//8

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__head,0)

            # Read data
            col = np.array(struct.unpack('d'*size, \
                                         f.read(self.__column))). \
                           reshape((self.__nz,self.__nn))

            # Close
            f.close()

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

################################################################################
################################################################################
################################################################################

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
          [None,'Get number of nodes in the x, y, and wavelength dimensions'], \
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
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract Stokes I at a particular column'], \
         'get_stokesq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular column'], \
         'get_stokesu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular column'], \
         'get_stokesv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular column'], \
         'get_linear_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a particular column'], \
         'get_stokes_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular column'], \
          'get_plane_stk': \
          [{'il': \
             'Coordinate in the wavelength dimension of the Stokes parameters to extract', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract Stokes parameters at a given wavelength position for ' + \
           'the whole field of view'], \
         'get_stokesi_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract Stokes I at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_stokesq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_stokesu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_stokesv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_linear_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokes_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular wavelength index for ' + \
           'the whole field of view']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D emergence file head
        '''
        debug = False
        try:
            # Get actual header
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

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_lambda,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                          f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
        except:
            raise

    def __get_gen_column(self,ix,iy,minl=None,maxl=None,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters column
        '''

        # Get size to read
        siz = self.__nl
        bsiz = siz*4

        # Output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Try geeting data
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
                        out[j] = np.array(struct.unpack('f'*siz, \
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

        # Failed
        except struct.error:

            # If the file is complete, the error is more severe,
            # let it crash
            if self.__complete:
                raise

            # Incomplete file, may be missing data
            else:

                # Warn
                msg = 'Could not read, may be due to the file being ' + \
                      'not complete'
                _error(msg,0)

                # Generate zeros
                for j in indx:
                    out[j] = np.zeros((self.__nl))

        # Others
        except:
            raise

        # Close file
        f.close()

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

    def _get_stokesq_column(self,ix,iy,minl=None,maxl=None,fractional=False):
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

        return self.__get_gen_column(ix,iy,minl,maxl,fractional,[1])[1]

    def _get_stokesu_column(self,ix,iy,minl=None,maxl=None,fractional=False):
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

        return self.__get_gen_column(ix,iy,minl,maxl,fractional,[2])[2]

    def _get_stokesv_column(self,ix,iy,minl=None,maxl=None,fractional=False):
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

        return self.__get_gen_column(ix,iy,minl,maxl,fractional,[3])[3]

    def _get_linear_column(self,ix,iy,minl=None,maxl=None,fractional=False):
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

    def _get_stokes_column(self,ix,iy,minl=None,maxl=None,fractional=False):
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

        iquv = self.__get_gen_column(ix,iy,minl,maxl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))


    def __get_gen_plane(self,il,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters plane
        '''

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

                # Try geeting data
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
                                out[j][ix,iy] = struct.unpack('f',f.read(4))[0]
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

                    # If the file is complete, the error is more severe, let it crash
                    if self.__complete:
                        raise

                    # Incomplete file, may be missing data
                    else:

                        # Warn
                        msg = 'Could not read, may be due to the file being not complete'
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



        # Units
        if fractional:
            if 0 in indx:
                out[0] *= self.__unit_trans
        else:
            for j in indx:
                out[j] *= self.__unit_trans

        # Close file
        f.close()

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

################################################################################
################################################################################
################################################################################

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
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, height, and wavelength dimensions'], \
         'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_ctri_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get intensity contribution function [SI] at a given column'], \
          'get_ctrq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get Stokes Q contribution function [SI] at a given column'], \
          'get_ctru_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get Stokes U contribution function [SI] at a given column'], \
          'get_ctrv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get Stokes V contribution function [SI] at a given column'], \
          'get_ctr_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary in wavelength [nm]', \
            'maxl': \
             'Upper boundary in wavelength [nm]'}, \
           'Get full Stokes vector contribution function [SI] at a given column']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D contribution function file head
        '''
        debug = False
        try:
            # Get actual header
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
        ''' Get number of positions in x, y, height, and wavelength axes
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

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_lambda,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                          f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
        except:
            raise

    def __get_gen_column(self,ix,iy,minl=None,maxl=None,indx=[0]):
        ''' Generic read of contribution function column
        '''

        # Get size to read
        siz = self.__nl*self.__nz
        bsiz = siz*4

        # Output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Try geeting data
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
                    out[j] = (np.array(struct.unpack('f'*siz, \
                                                    f.read(bsiz))). \
                                       reshape((self.__nz,self.__nl)))[:,::-1]
                else:

                    # Skip
                    f.seek(bsiz,1)

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

        # Failed
        except struct.error:

            # If the file is complete, the error is more severe,
            # let it crash
            if self.__complete:
                raise

            # Incomplete file, may be missing data
            else:

                # Warn
                msg = 'Could not read, may be due to the file being ' + \
                      'not complete'
                _error(msg,0)

                # Generate zeros
                for j in indx:
                    out[j] = np.zeros((self.__nz,self.__nl))

        # Others
        except:
            raise

        # Close file
        f.close()

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
        ''' Get full Stokes vector contribution function at a given column
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

################################################################################
################################################################################
################################################################################

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

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

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
          [None,'Get number of nodes in the x, y, and wavelength dimensions'], \
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
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
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
           'in optical depth at a particular wavelength index for ' + \
           'the whole field of view']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 1.5D height where optical depth is one file head
        '''
        debug = False
        try:
            # Get actual header
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
           #print(f'Head size {self.__head}')
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
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_lambda,0)
            omg = np.array(struct.unpack('d'*self.__nl, \
                                          f.read(8*self.__nl)))
            lam = 1e2/omg[::-1]
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
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

        # Get size to read
        siz = self.__nl
        bsiz = siz*4

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Try geeting data
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
                msg = 'Could not read, may be due to the file being ' + \
                      'not complete'
                _error(msg,0)

                # Generate zeros
                tau = np.zeros((self.__nl))

        # Others
        except:
            raise

        # Adjust wavelength
        if minl is not None:
            i = np.argmin(np.absolute(lam - minl))
            lam = lam[i:]
            tau = tau[i:]
        if maxl is not None:
            i = np.argmin(np.absolute(lam - maxl))
            lam = lam[:i+1]
            tau = tau[:i+1]

        # Close file
        f.close()

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

                # Try geeting data
                try:

                    # Get Stokes
                    if left > 0: f.seek(left,1)
                    tau[ix,iy] = struct.unpack('f',f.read(4))[0]
                    if right > 0: f.seek(right,1)

                # Reading error
                except struct.error:

                    # If the file is complete, the error is more severe, let it crash
                    if self.__complete:
                        raise

                    # Incomplete file, may be missing data
                    else:

                        # Warn
                        msg = 'Could not read, may be due to the file being not complete'
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

################################################################################
################################################################################
################################################################################

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
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_vars': \
          [None,'Get list of variables with available node results'], \
         'get_vars_units': \
          [None,'Get list of variables with available node results with their ' + \
                'corresponding units'], \
         'get_vars': \
          [None,'Get list of variables in the model atmosphere'], \
         'get_vars_units': \
          [None,'Get list of variables in the model atmosphere with ' + \
                'their corresponding units'], \
         'get_vars_alias': \
          [None,'Get list of variables in the model atmosphere with ' + \
                'their corresponding alias'], \
          'get_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minh': \
             'Lower boundary for output height', \
            'maxh': \
             'Upper boundary for output height', \
            'mint': \
             'Lower boundary for output optical depth', \
            'maxt': \
             'Upper boundary for output optical depth', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars_alias()}'}, \
           'Extract the model atmosphere for a particular column'], \
          'get_plane': \
          [{'iz': \
             'Coordinate in the height dimension of the atmospheric parameters to extract', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars_alias()}'}, \
           'Extract the model atmosphere for a particular height index for ' + \
           'the whole field of view']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert 3D model atmosphere
        '''
        try:

            # Get actual header
            f = open(self.__filename,'rb')
            f.seek(4,0)

            # Read precision
            size = int(struct.unpack('i',f.read(4))[0])
            if size == 4:
                self.__fmt = 'f'
                self.__byt = 4
            elif size == 8:
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
            self.__vars = ['h',r'$\tau_{\rm c}$',r'$\chi_{\rm c}$', \
                           'T',r'P$_{\rm g}$',r'$\rho$', \
                           r'B$_{\rm x}$',r'B$_{\rm y}$',r'B$_{\rm z}$', \
                           r'v$_{\rm x}$',r'v$_{\rm y}$',r'v$_{\rm z}$', \
                           r'v$_{\rm mi}$',r'P$_{\rm e}$', \
                           r'N$_{\rm e}$',r'N$_{\rm H}$', \
                           r'N$_{\rm H_{\rm a}}$', \
                           r'N$_{\rm H^-}$',r'N$_{\rm H_0}$',
                           r'N$_{\rm H_1}$',r'N$_{\rm H_2}$', \
                           r'N$_{\rm H_3}$',r'N$_{\rm H_4}$', \
                           r'N$_{\rm p^+}$']
            self.__alias = ['h','tau','chic','T','Pg','rho', \
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
        return self.__get_nxyz()

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

        # Size column
        size = self.__nz*self.__nvar

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__head + iy*self.__column + \
                       self.__ny*ix*self.__column,0)

            # Read data
            col = np.array(struct.unpack(self.__fmt*size, \
                                         f.read(self.__column))). \
                           reshape((self.__nvar,self.__nz))

            # Close
            f.close()

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

        except struct.error:

            # Not full data
            raise

        except:
            raise

        # Return column
        out = {}
        for i,v in enumerate(self.__alias):
            if v in ivar:
                out[v] = col[i,:]
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

        # Before and after
        left = iz*self.__byt
        right = (self.__nz - iz - 1)*self.__byt
        full = self.__nz*self.__byt

        # Prepare output
        out = {}
        for jvar in ivar:
            out[jvar] = np.empty((self.__nx,self.__ny))

        # Open file
        f = open(self.__filename,'rb')

        # Seek first data points for this block
        f.seek(self.__head,0)

        # Run over columns
        for ix in range(self.__nx):
            for iy in range(self.__ny):

                # Try geeting data
                try:

                    # For each variable in the model atmosphere
                    for ibar,bar in enumerate(self.__alias):

                        # Output
                        if bar in out:

                            # Skip left
                            if left > 0: f.seek(left,1)

                            # Read
                            out[bar][ix,iy] = struct.unpack(self.__fmt, \
                                                            f.read(self.__byt))[0]
                            # Skip right
                            if right > 0: f.seek(right,1)

                        else:

                            # Skip
                            f.seek(full,1)

                except struct.error:
                    raise
                except:
                    raise
        # Close
        f.close()

        return out

################################################################################
################################################################################
################################################################################

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
          [None,'Get number of nodes in the x, y, and wavelength dimensions'], \
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
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the intensity at a particular column'], \
         'get_stokesq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes Q at a particular column'], \
         'get_stokesu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes U at a particular column'], \
         'get_stokesv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract Stokes V at a particular column'], \
         'get_linear_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'fractional':
             'If normalized to the intensity'}, \
           'Extract linear polarization at a particular column'], \
         'get_stokes_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
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
           'Extract Stokes I at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_stokesq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes Q at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_stokesu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes U at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_stokesv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract Stokes V at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_linear_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract total linear polarization at a particular wavelength ' + \
           'index for the whole field of view'], \
         'get_stokes_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract', \
            'fractional': \
             'True to normalize to intensity, [SI] otherwise]'}, \
           'Extract the full Stokes vector at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_sigma': \
          [None,'Get the sigma for Stokes parameters'], \
         'get_sigmai_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the intensity sigma at a particular column'], \
         'get_sigmaq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes Q sigma at a particular column'], \
         'get_sigmau_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes U sigma at a particular column'], \
         'get_sigmav_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes V sigma at a particular column'], \
         'get_sigma_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes sigma at a particular column'], \
         'get_sigmai_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the intensity sigma at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_sigmaq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes Q sigma at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_sigmau_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes U sigma at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_sigmav_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes V sigma at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_sigma_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes sigma at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diff': \
          [None,'Get the diffuse light Stokes parameters'], \
         'get_diffi_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the intensity diffuse light at a particular column'], \
         'get_diffq_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes Q diffuse light at a particular column'], \
         'get_diffu_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes U diffuse light at a particular column'], \
         'get_diffv_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes V diffuse light at a particular column'], \
         'get_diff_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Extract the Stokes diffuse light at a particular column'], \
         'get_diffi_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the intensity diffuse light at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diffq_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes Q diffuse light at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diffu_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes U diffuse light at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diffv_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes V diffuse light at a particular wavelength index for ' + \
           'the whole field of view'], \
         'get_diff_plane': \
          [{'il': \
             'Coordinate in the wavelength dimension to extract'}, \
           'Extract the Stokes diffuse light at a particular wavelength index for ' + \
           'the whole field of view']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert inversion output file head
        '''
        try:

            # Get actual header
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
        self.__jump_to_lambda = 4*8
        self.__head = self.__jump_to_lambda + self.__nl*8
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
            self.__c_to_los = 0
            self.__c_los = 0
        # Variable LOS
        elif self.__info[1] == 1:
            self.__to_los = 0
            self.__c_to_los = 0
            self.__c_los = 16

        # Profiles
        self.__c_to_stk = self.__c_to_los + self.__c_los
        # If intensity
        if self.__info[0] == 0:
            self.__c_stk = self.__nl*8
        # If polarization
        elif self.__info[0] == 1:
            self.__c_stk = self.__nl*8*4

        # No sigma
        if self.__info[2] == 0:
            self.__to_sigma = self.__to_data
            self.__c_to_sigma = 0
            self.__c_sigma = 0
        # If constant sigma intensity
        elif self.__info[2] == 1 and self.__info[0] == 0:
            self.__to_sigma = self.__to_data
            self.__to_data += 8
            self.__c_to_sigma = 0
            self.__c_sigma = 0
        # If constant sigma polarization
        elif self.__info[2] == 1 and self.__info[0] == 1:
            self.__to_sigma = self.__to_data
            self.__to_data += 8*4
            self.__c_to_sigma = 0
            self.__c_sigma = 0
        # If variable sigma intensity
        elif self.__info[2] == 2 and self.__info[0] == 0:
            self.__to_sigma = self.__to_data
            self.__to_data += 8*self.__nl
            self.__c_to_sigma = 0
            self.__c_sigma = 0
        # If variable sigma polarization
        elif self.__info[2] == 2 and self.__info[0] == 1:
            self.__to_sigma = self.__to_data
            self.__to_data += 8*self.__nl*4
            self.__c_to_sigma = 0
            self.__c_sigma = 0
        # If pixel constant sigma intensity
        elif self.__info[2] == 3 and self.__info[0] == 0:
            self.__to_sigma = 0
            self.__c_to_sigma = self.__c_to_stk + self.__c_stk
            self.__c_sigma = 8
        # If pixel constant sigma polarization
        elif self.__info[2] == 3 and self.__info[0] == 1:
            self.__to_sigma = 0
            self.__c_to_sigma = self.__c_to_stk + self.__c_stk
            self.__c_sigma = 8*4
        # If pixel variable sigma intensity
        elif self.__info[2] == 4 and self.__info[0] == 0:
            self.__to_sigma = 0
            self.__c_to_sigma = self.__c_to_stk + self.__c_stk
            self.__c_sigma = 8*self.__nl
        # If pixel variable sigma polarization
        elif self.__info[2] == 4 and self.__info[0] == 1:
            self.__to_sigma = 0
            self.__c_to_sigma = self.__c_to_stk + self.__c_stk
            self.__c_sigma = 8*self.__nl*4

        # No diffuse light
        if self.__info[3] == 0:
            self.__to_diff = self.__to_data
            self.__c_to_diff = 0
            self.__c_diff = 0
        # If constant diffuse light intensity
        elif self.__info[3] == 1:
            self.__to_diff = self.__to_data
            self.__to_data += 8*self.__nl
            self.__c_to_diff = 0
            self.__c_diff = 0
        # If constant diffuse light polarization
        elif self.__info[3] == 2:
            self.__to_diff = self.__to_data
            self.__to_data += 8*4
            self.__c_to_diff = 0
            self.__c_diff = 0
        # If pixel diffuse light intensity
        elif self.__info[3] == 3:
            self.__to_diff = 0
            self.__c_to_diff = self.__c_to_stk + self.__c_stk + \
                               self.__c_sigma
            self.__c_diff = self.__nl*8
        # If pixel diffuse light polarization
        elif self.__info[3] == 4:
            self.__to_diff = 0
            self.__c_to_diff = self.__c_to_stk + self.__c_stk + \
                               self.__c_sigma
            self.__c_diff = self.__nl*8*4

        # Size of a column
        self.__column = self.__c_los + self.__c_stk + self.__c_sigma + self.__c_diff
        self.__c_from_los = self.__c_stk + self.__c_sigma + self.__c_diff
        self.__c_from_stk = self.__c_sigma + self.__c_diff
        self.__c_from_sigma = self.__c_diff
        self.__c_from_diff = 0

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
            _error('Only intensity file',0)
        else:
            _error('Full Stokes file',0)
        if self.__info[1] == 0:
            _error('Single LOS',0)
        else:
            _error('Pixelwise LOS',0)
        if self.__info[2] == 0:
            _error('Constant sigma',0)
        elif self.__info[2] == 1:
            _error('Wavelength dependent constant sigma',0)
        elif self.__info[2] == 2:
            _error('Constant pixelwise sigma',0)
        elif self.__info[2] == 3:
            _error('Wavelength dependent pixelwise sigma',0)
        if self.__info[3] == 0:
            _error('Constant only intensity diffuse light profile',0)
        elif self.__info[3] == 1:
            _error('Constant full Stokes diffuse light profile',0)
        elif self.__info[3] == 2:
            _error('Pixelwise only intensity diffuse light profile',0)
        elif self.__info[3] == 3:
            _error('Pixelwise full Stokes diffuse light profile',0)

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
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_lambda,0)
            lam = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
        except:
            raise

    def __get_los_ct(self):
        ''' Get constant los
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_los,0)
            th = struct.unpack('d',f.read(8))[0]
            ph = struct.unpack('d',f.read(8))[0]
            f.close()
        except struct.error:
            raise
        except:
            raise
        return [th,ph]

    def __get_los_column(self,ix,iy):
        ''' Get los at a column
        '''
        try:
            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__jump_to_data + iy*self.__column + \
                               self.__ny*ix*self.__column,0)

            # Jump to los
            if self.__c_to_los > 0: f.seek(self.__c_to_los,1)

            # Get LOS
            th = struct.unpack('d',f.read(8))[0]
            ph = struct.unpack('d',f.read(8))[0]

            f.close()

        except struct.error:
            raise
        except:
            raise
        return [th,ph]

    def __get_los_plane(self):
        ''' Get los for the plane
        '''

        # Output
        los = np.empty((self.__nx,self.__ny,2))

        try:
            # Open file
            f = open(self.__filename,'rb')

            # Seek first data point
            f.seek(self.__jump_to_data,0)

            # For each column
            for ix in range(self.__nx):
                for iy in range(self.__ny):

                    # Jump to los
                    if self.__c_to_los > 0: f.seek(self.__c_to_los,1)

                    # Get LOS
                    los[ix,iy,:] = np.array(struct.unpack('dd',f.read(16)))

                    # End column
                    if self.__c_from_los > 0: f.seek(self.__c_from_los,1)

            f.close()

        except struct.error:
            raise
        except:
            raise

        return los

    def _get_los(self):
        ''' Get LOS if full constant
        '''
        # Constant
        if self.__info[1] == 0:
            return np.array(self.__get_los_ct())
        else:
            return self.__get_los_plane()

    def _get_los_column(self,ix=None,iy=None):
        ''' Get LOS at a given column
        '''
        # Constant
        if self.__info[1] == 0:
            return np.array(self.__get_los_ct())
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
            return np.array(self.__get_los_column(ix,iy))

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
            return __get_los_plane()

    def __get_gen_column(self,ix,iy,ilvar,irvar,pol,minl=None,maxl=None,nl=None,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters column
        '''

        # Output
        out = [None,None,None,None]

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Sizes
        siz = nl
        bsiz = siz*8

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__to_data + iy*self.__column + \
                          self.__ny*ix*self.__column,0)

            # Skip to the left?
            if ilvar > 0: f.seek(ilvar,1)

            # Intensity
            if 0 in indx or fractional:

                # Get intensity
                stkI = np.array(struct.unpack('d'*siz, \
                                              f.read(bsiz)))

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
                if pol:

                    # To output
                    if j in indx:

                        # Read Stokes
                        out[j] = np.array(struct.unpack('d'*siz, \
                                                         f.read(bsiz)))
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

            # Skip to the right?
            if irvar > 0: f.seek(irvar,1)

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

        # Failed
        except struct.error:
            raise
        # Others
        except:
            raise

        # Close file
        f.close()

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

        return self.__get_gen_column(ix,iy,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,minl,maxl,self.__nl,False,[0])[0]

    def _get_stokesq_column(self,ix,iy,minl=None,maxl=None,fractional=False):
        ''' Get Stokes Q profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
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

        return self.__get_gen_column(ix,iy,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,minl,maxl,self.__nl,fractional,[1])[1]

    def _get_stokesu_column(self,ix,iy,minl=None,maxl=None,fractional=False):
        ''' Get Stokes U profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
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

        return self.__get_gen_column(ix,iy,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,minl,maxl,self.__nl,fractional,[2])[2]

    def _get_stokesv_column(self,ix,iy,minl=None,maxl=None,fractional=False):
        ''' Get Stokes V profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
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

        return self.__get_gen_column(ix,iy,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,minl,maxl,self.__nl,fractional,[3])[3]

    def _get_linear_column(self,ix,iy,minl=None,maxl=None,fractional=False):
        ''' Get Stokes linear polarization profile at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
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

        qu = self.__get_gen_column(ix,iy,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,minl,maxl,self.__nl,fractional,[1,2])
        return np.sqrt(qu[1]*qu[1] + qu[2]*qu[2])

    def _get_stokes_column(self,ix,iy,minl=None,maxl=None,fractional=False):
        ''' Get Stokes parameter at a given column
        '''

        # Mode?
        if self.__info[0] == 0:
           _error('The file is only intensity',1)
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

        iquv = self.__get_gen_column(ix,iy,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,minl,maxl,self.__nl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))


    def __get_gen_plane(self,il,ilvar,irvar,pol,nl,fractional=False,indx=[0]):
        ''' Generic read of Stokes parameters plane
        '''

        # Get size to read
        left = il*8
        right = (nl - il - 1)*8
        full = nl*8
        abort = False

        # Output
        out = [None,None,None,None]

        # For each index requested
        for j in indx:
            out[j] = np.empty((self.__nx,self.__ny))

        # Open file
        f = open(self.__filename,'rb')

        # Seek to data
        f.seek(self.__to_data,0)

        # For each column
        for ix in range(self.__nx):
            for iy in range(self.__ny):

                # Try geeting data
                try:

                    # Left vars
                    if ilvar > 0: f.seek(ilvar,1)

                    # Intensity
                    if 0 in indx or fractional:

                        # Get intensity
                        if left > 0: f.seek(left,1)
                        stkI = struct.unpack('d',f.read(8))[0]
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
                        if pol:

                            # To output
                            if j in indx:

                                # Get Stokes
                                siz = 0
                                if left > 0: f.seek(left,1)
                                out[j][ix,iy] = struct.unpack('d',f.read(8))[0]
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

                    # Right vars
                    if irvar > 0: f.seek(irvar,1)

                # Reading error
                except struct.error:
                    raise
                except:
                    raise

        # Close file
        f.close()

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

        return self.__get_gen_plane(il,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,self.__nl,False,[0])[0]

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

        return self.__get_gen_plane(il,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,self.__nl,fractional,[1])[1]

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

        return self.__get_gen_plane(il,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,self.__nl,fractional,[2])[2]

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

        return self.__get_gen_plane(il,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,self.__nl,fractional,[3])[3]

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

        qu = self.__get_gen_plane(il,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,self.__nl,fractional,[1,2])
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

        iquv = self.__get_gen_plane(il,self.__c_to_stk,self.__c_from_stk,self.__info[0]==1,self.__nl,fractional,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def __no_sigma(self):
        ''' Message that there is no sigma
        '''
        __error('There is no sigma in this file',1)
        return None

    def __get_sigma_ct(self):
        ''' Get constant sigma
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_sigma,0)
            # Full constant
            if self.__info[2] == 1:
                if self.__info[0] == 0:
                    sig = [0.,0.,0.,0.]
                    sig[0] = struct.unpack('d',f.read(8))[0]
                else:
                    sig = struct.unpack('dddd',f.read(32))
            # Wavelength dependent constant
            elif self.__info[2] == 2:
                if self.__info[0] == 0:
                    sig = ['', \
                           np.zeros((self.__nl)), \
                           np.zeros((self.__nl)), \
                           np.zeros((self.__nl))]
                    sig[0] = np.array(struct.unpack('d'*self.__nl,f.read(8*self.__nl)))
                else:
                    sig = []
                    for i in range(4):
                        sig.append(np.array(struct.unpack('d'*self.__nl,f.read(8*self.__nl))))
            f.close()
        except struct.error:
            raise
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
            iquv = self.__get_sigma_plane(0,self.__c_to_sigma, \
                    self.__c_from_sigma,self.__info[0]==1,1,False,[0,1,2,3,4])
            return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))
        # Other
        else:
            _error('Sigma is pixelwise and non-constant, get a column with get_sigma_column() or a plane with get_sigma_plane()',1)
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

        return self.__get_gen_column(ix,iy,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,iminl,imaxl,nl,False,[0])[0]

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

        return self.__get_gen_column(ix,iy,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,iminl,imaxl,nl,False,[1])[1]

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

        return self.__get_gen_column(ix,iy,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,iminl,imaxl,nl,False,[2])[2]

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

        return self.__get_gen_column(ix,iy,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,iminl,imaxl,nl,False,[3])[3]

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

        iquv = self.__get_gen_column(ix,iy,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,iminl,imaxl,nl,False,[0,1,2,3])
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

        return self.__get_gen_plane(il,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,nl,False,[0])[0]

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

        return self.__get_gen_plane(il,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,nl,False,[1])[1]

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

        return self.__get_gen_plane(il,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,nl,False,[2])[2]

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

        return self.__get_gen_plane(il,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,nl,False,[3])[3]

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

        iquv = self.__get_gen_plane(il,self.__c_to_sigma,self.__c_from_sigma,self.__info[0]==1,nl,False,[0,1,2,3])
        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

    def __no_diff(self):
        ''' Message that there is no diffuse light
        '''
        __error('There is no diffuse light in this file',1)
        return None

    def __get_diff_ct(self):
        ''' Get constant diffuse light
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_diff,0)
            # Intensity
            if self.__info[3] == 1:
                diff = [np.array(struct.unpack('d'*self.__nl,f.read(8*self.__nl)))]
                diff.append(np.zeros((self.__nl)))
                diff.append(np.zeros((self.__nl)))
                diff.append(np.zeros((self.__nl)))
            # Polarized
            elif self.__info[3] == 2:
                diff = []
                for i in range(4):
                    diff.append(np.array(struct.unpack('d'*self.__nl,f.read(8*self.__nl))))
            f.close()
        except struct.error:
            raise
        except:
            raise
        return diff

    def _get_diff(self):
        ''' Get diff if full constant
        '''
        # No diff
        if self.__info[3] == 0: return self.__no_diff()
        # Constant
        if self.__info[3] == 0 or self.__info[3] == 1:
            return self.__get_diff_ct()
        # Other
        else:
            _error('Diffuse light is pixelwise, get a column with get_diff_column() or a plane with get_diff_plane()',1)
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

        return self.__get_gen_column(ix,iy,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,minl,maxl,self.__nl,False,[0])[0]

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

        return self.__get_gen_column(ix,iy,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,minl,maxl,self.__nl,False,[1])[1]

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

        return self.__get_gen_column(ix,iy,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,minl,maxl,self.__nl,False,[2])[2]

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

        return self.__get_gen_column(ix,iy,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,minl,maxl,self.__nl,False,[3])[3]

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

        iquv = self.__get_gen_column(ix,iy,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,minl,maxl,self.__nl,False,[0,1,2,3])
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

        return self.__get_gen_plane(il,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,self.__nl,False,[0])[0]

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

        return self.__get_gen_plane(il,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,self.__nl,False,[1])[1]

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

        return self.__get_gen_plane(il,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,self.__nl,False,[2])[2]

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

        return self.__get_gen_plane(il,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,self.__nl,False,[3])[3]

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

        iquv = self.__get_gen_plane(il,self.__c_to_diff,self.__c_from_diff,self.__info[3]==4,self.__nl,False,[0,1,2,3])

        return np.stack((iquv[0],iquv[1],iquv[2],iquv[3]))

################################################################################
################################################################################
################################################################################

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

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

        # Method
        self.__methods = { \
         'get_filename': \
          [None,'Get name of the read file'], \
         'get_polarized': \
          [None,'Get if the inversion included polarization'], \
         'get_nx': \
          [None,'Get number of nodes in the x dimension'], \
         'get_ny': \
          [None,'Get number of nodes in the y dimension'], \
         'get_nz': \
          [None,'Get number of heights in the model atmosphere'], \
         'get_nxy': \
          [None,'Get number of nodes in the x and y dimensions'], \
         'get_nxyz': \
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, height, and wavelength dimensions'], \
         'get_vars': \
          [None,'Get list of variables with available node results'], \
         'get_vars_units': \
          [None,'Get list of variables with available node results with their ' + \
                'corresponding units'], \
         'get_vars_atmo': \
          [None,'Get list of variables in the model atmosphere'], \
         'get_vars_atmo_units': \
          [None,'Get list of variables in the model atmosphere with ' + \
                'their corresponding units'], \
          'get_lambda': \
          [{'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]'}, \
           'Get wavelengths in [nm]'], \
          'get_column': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars()}'}, \
           'Extract the result of the inversion at a particular column'], \
          'get_column_atmo': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minh': \
             'Lower boundary for output optical depth', \
            'maxh': \
             'Upper boundary for output optical depth', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars_atmo()}'}, \
           'Extract the model atmosphere for a particular column'], \
          'get_column_rf': \
          [{'ix': \
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract the response function at a particular column'], \
         'get_plane_chi': \
          [None,'Get the value of the initial and final merit function for ' + \
                'the whole field of view'], \
          'get_plane_stk': \
          [{'il': \
             'Coordinate in the wavelength dimension of the Stokes parameters to extract', \
            'var': \
             'List of variables to include in the output'}, \
           'Extract Stokes parameters at a given wavelength position for ' + \
           'the whole field of view'], \
          'get_plane_atmo': \
          [{'iz': \
             'Coordinate in the height dimension of the atmospheric parameters to extract', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars_atmo()}'}, \
           'Extract the model atmosphere for a particular height index for ' + \
           'the whole field of view'], \
          'get_node': \
          [{'var': \
             'Variables for which to extract the node information'}, \
           'Extract the node information (full cube) for a given variable']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert inversion output file head
        '''
        debug = False
        try:
            # Get actual header
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
              self.__vars_atmo = ['tau','ltau','T','Pg', \
                                  'Bx','By','Bz','vx', \
                                  'vy','vz','vturb','ne', \
                                  'nHT','nHa','nH(0)','nH(1)', \
                                  'nH(2)','nH(3)','nH(4)','np', \
                                  'J10','Re{J11}','Im{J11}', \
                                  'J20','Re{J21}','Im{J21}', \
                                  'Re{J22}','Im{J22}','f']
              self.__vars_atmo_units = ['','','K','dyn cm^-2','G', \
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
              self.__vars_atmo = ['tau','ltau','T','Pg', \
                                  'Bx','By','Bz','vx', \
                                  'vy','vz','vturb','ne', \
                                  'nHT','nHa','nH(0)','nH(1)', \
                                  'nH(2)','nH(3)','nH(4)','np','f']
              self.__vars_atmo_units = ['','','K','dyn cm^-2','G', \
                                        'G','G','km s^-1','km s^-1', \
                                        'km s^-1','km s^-1','cm^-3', \
                                        'cm^-3','cm^-3','cm^-3', \
                                        'cm^-3','cm^-3','cm^-3', \
                                        'cm^-3','cm^-3','']

            # Atmosphere pixel size                          v f_diff
            self.__s_atmo_c = self.__nz*self.__nvar_atmo*4 + 4

            # Atmosphere size
            __s_atmo = self.__nx*self.__ny*self.__s_atmo_c

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
                self.__vars_add_units = ['','','J m^-2 s^-1 sr Hz^-1', \
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
                self.__vars_add_units = ['','','J m^-2 s^-1 sr Hz^-1', \
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
                self.__inv_flag.append(struct.unpack('i', f.read(4))[0] > 0)

                # Number of nodes
                self.__inv_node.append(struct.unpack('i', f.read(4))[0])

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
                    self.__inv_vnode.append(struct.unpack('i',f.read(4))[0])
                    if debug: print('Get nodes',self.__inv_vnode[-1])

                    # Add sizes
                    __s_rf_head += 8

                    # If varying nodes
                    if self.__inv_vnode[-1] > 0:
                        if self.__polarization:
                            self.__s_rf_c += (4 + 16*self.__nl)*self.__inv_vnode[-1]
                        else:
                            self.__s_rf_c += (4 + 4*self.__nl)*self.__inv_vnode[-1]

                # Jump to RF
                self.__jump_to_rf = self.__jump_to_res + __s_res + __s_rf_head
                self.__is_RF = True

            except struct.error:

                self.__is_RF = False

            except:

                raise

        else:

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
        ''' Get number of positions in x, y, height, and wavelength axes
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
        ''' Get variables available in the atmospheric model with units
        '''
        out = []
        for var,uni in zip(self.__vars_atmo,self.__vars_atmo_units):
            out.append(var+' ['+uni+']')
        return out

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_lambda,0)
            lam = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
            return lam
        except struct.error:
            raise
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

        # Get column size
        bsiz = self.__s_res_c
        siz = bsiz//4

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__jump_to_res + iy*bsiz + self.__ny*ix*bsiz,0)

            # Chi_ori
            chi0 = struct.unpack('f',f.read(4))[0]
            chi = struct.unpack('f',f.read(4))[0]

            # Stokes
            if self.__polarization:
                stokes_ob = np.array(struct.unpack('f'*self.__nl*4, \
                                                   f.read(16*self.__nl))). \
                                    reshape((4,self.__nl))
                stokes_fi = np.array(struct.unpack('f'*self.__nl*4, \
                                                   f.read(16*self.__nl))). \
                                    reshape((4,self.__nl))
            else:
                stokes_ob = np.array(struct.unpack('f'*self.__nl, \
                                                   f.read(4*self.__nl))). \
                                    reshape((1,self.__nl))
                stokes_fi = np.array(struct.unpack('f'*self.__nl, \
                                                   f.read(4*self.__nl))). \
                                    reshape((1,self.__nl))

            # Adjust wavelength
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                stokes_ob = stokes_ob[:,i:]
                stokes_fi = stokes_ob[:,i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                stokes_ob = stokes_ob[:,:i+1]
                stokes_fi = stokes_ob[:,:i+1]

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
                nodes[var][0:2,:] = np.array(struct.unpack('f'*node*2, \
                                                           f.read(8*node))). \
                                            reshape((2,node))
                # Error?
                if self.__inv_flag[ipar]:
                    nodes[var][2,:] = np.array(struct.unpack('f'*node, \
                                                             f.read(4*node)))
                else:
                    nodes[var][2,:] = 0e0


            # Close
            f.close()

        # Failed
        except struct.error:

            # If the file is complete, the error is more severe, let it crash
            if self.__complete:
                raise

            # Incomplete file, may be missing data
            else:

                # Close and warn
                f.close()
                msg = 'Could not read, may be due to the file being not complete'
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

        # Get column size
        bsiz = self.__s_atmo_c
        siz = bsiz//4

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__jump_to_atmo + iy*bsiz + self.__ny*ix*bsiz,0)

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

            # Adjust height
            if minh is not None:
                i = np.argmin(np.absolute(col[0,:] - minh))
                col = col[:,:i+1]
            if maxh is not None:
                i = np.argmin(np.absolute(col[0,:] - maxh))
                col = col[:,i:]

        except struct.error:

            # If the file is complete, the error is more severe, let it crash
            if self.__complete:
                raise

            # Incomplete file, may be missing data
            else:

                # Close and warn
                f.close()
                msg = 'Could not read, may be due to the file being not complete'
                _error(msg,0)

                # Generate zeros
                col = np.zeros((self.__nvar_atmo,self.__nz))
                f_diff = 0

        except:
            raise

        # Return column
        out = {}
        if 'tau' in ivar:
            out['tau'] = col[0,:]
        if 'ltau' in ivar:
            out['ltau'] = np.log10(col[0,:])
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

        # Get column size
        bsiz = self.__s_rf_c
        siz = bsiz//4

        # Need lambda?
        if minl is not None or maxl is not None:
            lam = self._get_lambda()

        # Output
        out = {}

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__jump_to_rf + iy*bsiz + self.__ny*ix*bsiz,0)

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
                red[jvar] = {'H': np.array(struct.unpack('f'*vnode, \
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
                        np.array(struct.unpack('f'*vnode*4*nl, \
                                               f.read(16*nl*vnode))). \
                                 reshape((vnode,4,self.__nl))
                else:
                    red[jvar]['RF'] = \
                        np.array(struct.unpack('f'*vnode*nl, \
                                               f.read(4*nl*vnode))). \
                                 reshape((vnode,1,self.__nl))

                # Add to output?
                if self.__vars_dic[jvar] in ivar:
                    out[self.__vars_dic[jvar]] = red[jvar]

            # Close file
            f.close()

        # Failed
        except struct.error:

            # If the file is complete, the error is more severe, let it crash
            if self.__complete:
                raise

            # Incomplete file, may be missing data
            else:

                # Close and warn
                f.close()
                msg = 'Could not read, may be due to the file being not complete'
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
                        red[jvar] = {'RF': np.zeros((vnode,4,self.__nl))}
                    else:
                        red[jvar] = {'RF': np.zeros((vnode,1,self.__nl))}

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

                    chi[ix,iy,:] = np.array(struct.unpack('ff',f.read(8)))

                except struct.error:
                    if self.__complete:
                        raise
                    else:
                        msg = 'Could not read, may be due to the file ' + \
                              'being not complete'
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
        ''' Get observed or fitted Stokes parameters for a given wavelength
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

        # Get column size
        bsiz = self.__s_res_c

        # Before and after
        before = 8
        left = il*4
        right = (self.__nl - il - 1)*4
        full = self.__nl*4
        if self.__polarization:
            after = bsiz - 8 - self.__nl*32
            stko = ['Io','Qo','Uo','Vo']
            stkf = ['If','Qf','Uf','Vf']
        else:
            after = bsiz - 8 - self.__nl*8
            stko = ['Io']
            stkf = ['If']

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

                # Try geeting data
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
                            out[stk][ix,iy] = struct.unpack('f',f.read(4))[0]
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
                            out[stk][ix,iy] = struct.unpack('f',f.read(4))[0]
                            # Skip right
                            if right > 0: f.seek(right,1)

                        else:

                            # Skip
                            f.seek(full,1)

                except struct.error:
                    if self.__complete:
                        raise
                    else:
                        msg = 'Could not read, may be due to the file ' + \
                              'being not complete'
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

                # Try geeting data
                try:

                    # For each variable in the model atmosphere
                    for ibar,bar in enumerate(self.__vars_atmo):

                        # Skip ltau
                        if bar == 'ltau': continue

                        # If tau
                        if bar == 'tau':

                            # If either tau or ltau in the list
                            if bar in out or 'ltau' in out:

                                # Left
                                if left > 0: f.seek(left,1)

                                # Get tau
                                aux = struct.unpack('f',f.read(4))[0]

                                # Skip right
                                if right > 0: f.seek(right,1)

                                # If output tau
                                if bar in out:
                                    out[bar] = aux
                                if 'ltau' in out:
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
                                out[bar] = struct.unpack('f',f.read(4))[0]

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
                              out[bar][ix,iy] = struct.unpack('f',f.read(4))[0]
                              # Skip right
                              if right > 0: f.seek(right,1)

                          else:

                              # Skip
                              f.seek(full,1)

                except struct.error:
                    if self.__complete:
                        raise
                    else:
                        msg = 'Could not read, may be due to the file ' + \
                              'being not complete'
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
        if evar not in self.__vars:
            _error('The requested variable ' + evar + \
                   ' is not available, ' + \
                   'check with get_vars',1)
            return None

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

                # Try geeting data
                try:

                    # Skip before and left
                    f.seek(before + left)

                    # Get nodes
                    for col in range(nsiz):
                        out[ix,iy,col,:] = \
                                 np.array(struct.unpack('f'*nnode, \
                                                        f.read(4*nnode)))
                except struct.error:
                    if self.__complete:
                        raise
                    else:
                        msg = 'Could not read, may be due to the file ' + \
                              'being not complete'
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

################################################################################
################################################################################
################################################################################

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

        #  Transformation to SI
        self.__unit_trans = 1e0/299792458e5

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
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_dims': \
          [None,' Get number of positions in x, y, and height axes and entries']}
        # Collisions
        if self.__cols:
            self.__methods['get_type'] = \
              [None,'Get the type of collisional rates in the file']
            self.__methods['get_nentry'] = \
              [None,'Get number of collisional entries']
            self.__methods['get_column'] = \
              [{'ix': \
                 'Coordinate in the x dimension of the column to extract', \
                'iy': \
                 'Coordinate in the y dimension of the column to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the collisional rates for a particular column']
            self.__methods['get_plane'] = \
              [{'iz': \
                 'Coordinate in the height dimension of the atmospheric parameters to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the collisional rates for a particular height index for ' + \
               'the whole field of view']
        # Damping
        else:
            self.__methods['get_nentry'] = \
              [None,'Get number of damping entries']
            self.__methods['get_column'] = \
              [{'ix': \
                 'Coordinate in the x dimension of the column to extract', \
                'iy': \
                 'Coordinate in the y dimension of the column to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the dampings for a particular column']
            self.__methods['get_plane'] = \
              [{'iz': \
                 'Coordinate in the height dimension of the atmospheric parameters to extract', \
                'ie': \
                 'List of entries to include in the output'}, \
               'Extract the dampings for a particular height index for ' + \
               'the whole field of view']

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert collisional/damping file head
        '''
        try:

            # Get actual header
            f = open(self.__filename,'rb')

            # Skip 3 first letter
            f.seek(3,0)

            # Read the fourth for the type of collisions
            self.__type = f.read(1).decode('utf-8')
            self.__cols = self.__type == 'l' or self.__type == 't'

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
        ''' Get number of positions in x, y, and height axes and entries
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

        # Size column
        size = self.__column//4

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__head + iy*self.__column + \
                       self.__ny*ix*self.__column,0)

            # Read data
            col = np.array(struct.unpack('f'*size, \
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
                msg = 'Could not read, may be due to the file being ' + \
                      'not complete'
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
            col = np.delete(col,np.array(todel,dtype='int32'),axis=0)

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

                # Try geeting data
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
                        out[ix,iy,ibar] = struct.unpack('f',f.read(4))[0]

                        # Skip right
                        if right > 0: f.seek(right,1)

                except struct.error:

                    # If the file is complete, the error is more severe, let it crash
                    if self.__complete:
                        raise

                    # Incomplete file, may be missing data
                    else:

                        # Warn
                        msg = 'Could not read, may be due to the file being not complete'
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
            out = np.delete(out,np.array(todel,dtype='int32'),axis=2)

        # If collisions, units factor
        if self.__cols: out *= 1e8

        return out

################################################################################
################################################################################
################################################################################

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
       #self.__unit_trans = 1e0/299792458e3
        #  Transformation to CGS
        self.__unit_trans = 1e0/299792458e2

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
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_nl': \
          [None,'Get number of wavelengths'], \
         'get_dims': \
          [None,'Get number of nodes in the x, y, height, and wavelength dimensions'], \
         'get_vars': \
          [None,'Get list of continuum variables'], \
         'get_vars_alias': \
          [None,'Get list of continuum variables with their aliases'], \
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
             'Coordinate in the x dimension of the column to extract', \
            'iy': \
             'Coordinate in the y dimension of the column to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars_alias()}'}, \
           'Extract the continuum variables at a particular column'], \
          'get_plane': \
          [{'iz': \
             'Coordinate in the z dimension of the plane to extract', \
            'minl': \
             'Lower boundary for output wavelength [nm]', \
            'maxl': \
             'Upper boundary for output wavelength [nm]', \
            'var': \
             'List of variables to include in the output (see the ' + \
             'available ones with get_vars_alias()}'}, \
           'Extract the continuum variables at a particular plane']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert background continuum file head
        '''
        debug = False
        try:
            # Get actual header
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
        ''' Get number of positions in x, y, height, and wavelength axes
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

    def _get_lambda(self,minl=None,maxl=None):
        ''' Get lambda from file
        '''
        try:
            f = open(self.__filename,'rb')
            f.seek(self.__jump_to_lambda,0)
            lam = np.array(struct.unpack('d'*self.__nl, \
                                         f.read(8*self.__nl)))
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
            f.close()
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
            lam = self._get_lambda()

        # Size
        siz = self.__column//4

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__head + iy*self.__column + \
                       self.__ny*ix*self.__column,0)

            # Read data
            col = np.array(struct.unpack('f'*siz, \
                                         f.read(self.__column))). \
                           reshape((self.__nz,self.__nvar,self.__nl))

            # Close
            f.close()

            # Adjust wavelength
            if minl is not None:
                i = np.argmin(np.absolute(lam - minl))
                lam = lam[i:]
                col = col[:,:,i:]
            if maxl is not None:
                i = np.argmin(np.absolute(lam - maxl))
                lam = lam[:i+1]
                col = col[:,:,:i+1]

        except struct.error:

            # If the file is complete, the error is more severe, let it crash
            if self.__complete:
                raise

            # Incomplete file, may be missing data
            else:

                # Close and warn
                f.close()
                msg = 'Could not read, may be due to the file being not complete'
                _error(msg,0)

                # Generate zeros
                col = np.zeros((self.__nz,self.__nvar,self.__nl))

        except:
            raise

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
            lam = self._get_lambda()

        # Size
        left = iz*self.__nvar*self.__nl*4
        right = (self.__nz - iz - 1)*self.__nvar*self.__nl*4
        siz = self.__nvar*self.__nl
        bsiz = siz*4
        abort = False

        # Prepare reading
        col = np.empty((self.__nx,self.__ny,self.__nvar,self.__nl))

        # Open file
        f = open(self.__filename,'rb')

        # Seek first data points for this block
        f.seek(self.__head,0)

        # Run over columns
        for ix in range(self.__nx):
            for iy in range(self.__ny):

                # Try geeting data
                try:

                    # Skip left
                    if left > 0: f.seek(left,1)

                    # Read
                    col[ix,iy,:,:] = np.array(struct.unpack('f'*siz, \
                                                            f.read(bsiz))). \
                                             reshape((self.__nvar,self.__nl))

                    # Skip right
                    if right > 0: f.seek(right,1)

                except struct.error:

                    # If the file is complete, the error is more severe, let it crash
                    if self.__complete:
                        raise

                    # Incomplete file, may be missing data
                    else:

                        # Warn
                        msg = 'Could not read, may be due to the file being not complete'
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

################################################################################
################################################################################
################################################################################

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
          [None,'Get number of nodes in the x, y, and height dimensions'], \
         'get_dims': \
          [None,' Get number of positions in x, y, and height axes and entries'], \
         'get_column': \
          [{'ix': \
            'Coordinate in the x dimension of the column to extract', \
            'iy': \
            'Coordinate in the y dimension of the column to extract', \
            'ie': \
            'List of entries to include in the output'}, \
            'Extract the populations/departures for a particular column'], \
          'get_plane': \
          [{'iz': \
            'Coordinate in the height dimension of the atmospheric parameters to extract', \
            'ie': \
            'List of entries to include in the output'}, \
            'Extract the populations/departures for a particular height index for ' + \
            'the whole field of view']}

    def _get_help(self):
        ''' Return methods dictionary
        '''
        return self.__methods

    def __head(self):
        ''' Reads hanlert collisional/damping file head
        '''
        try:

            # Get actual header
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
        ''' Get number of positions in x, y, and height axes and entries
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

        # Size column
        size = self.__column//4

        # Try geeting data
        try:

            # Open file
            f = open(self.__filename,'rb')

            # Seek first data points for this column
            f.seek(self.__head + iy*self.__column + \
                       self.__ny*ix*self.__column,0)

            # Read data
            col = np.array(struct.unpack('f'*size, \
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
                msg = 'Could not read, may be due to the file being ' + \
                      'not complete'
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
            col = np.delete(col,np.array(todel,dtype='int32'),axis=1)

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

                # Try geeting data
                try:

                    # Skip left
                    if left > 0: f.seek(left,1)

                    # Read
                    out[ix,iy,:] = np.array(struct.unpack('f'*self.__nn, \
                                                          f.read(4*self.__nn)))
                    # Skip right
                    if right > 0: f.seek(right,1)

                except struct.error:

                    # If the file is complete, the error is more severe, let it crash
                    if self.__complete:
                        raise

                    # Incomplete file, may be missing data
                    else:

                        # Warn
                        msg = 'Could not read, may be due to the file being not complete'
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
            out = np.delete(out,np.array(todel,dtype='int32'),axis=2)

        # Units
        if self.__type == 'p':
            out *= self.__unit_trans

        return out

################################################################################
################################################################################
################################################################################

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

    def __get_class(self):
        ''' Identify the type of file and link the suitable class
        '''

        # TODO TODO
        '''
  'CLEe': 'Emergent Stokes parameters from CLE synthesis', \
  'CLEt': 'Optical depth from CLE synthesis', \
  'MRC': 'Maximum relative change from 1.5D synthesis', \
  'sp': 'Solution file with polarization', \
  'si': 'Solution file without polarization', \
  'bo': 'Stokes parameters in the quadrature in 1D synthesis', \
  'ko': 'Frequency dependent radiation field tensors in 1D synthesis', \
  'ct': 'Term to term collisional rates from 1D synthesis', \
  'cl': 'Level to level collisional rates from 1D synthesis', \
  'da': 'Damping parameters from 1D synthesis', \
  'ba': 'Background continuum from 1D synthesis'}
        '''
        # TODO TODO

        # Possible labels
        labels = {'invo': 'Inversion Result file', \
                  'invi': 'Inversion input file', \
                  '2Dbe': 'Emergent Stokes parameters from 1.5D synthesis', \
                  '2Dbc': 'Contribution function from 1.5D synthesis', \
                  '2Dbt': 'Height for optical depth unity from 1.5D synthesis', \
                  '2Dct': 'Term to term collisional rates from 1.5D synthesis', \
                  '2Dcl': 'Level to level collisional rates from 1.5D synthesis', \
                  '2Dda': 'Damping parameters from 1.5D synthesis', \
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
                  'bo': 'Stokes parameters in the quadrature in 1D synthesis', \
                  'ko': 'Frequency dependent radiation field tensors in 1D synthesis', \
                  'be': 'Emergent Stokes parameters from 1D synthesis', \
                  'bc': 'Contribution function from 1D synthesis', \
                  'bt': 'Height for optical depth unity from 1D synthesis', \
                  'ct': 'Term to term collisional rates from 1D synthesis', \
                  'cl': 'Level to level collisional rates from 1D synthesis', \
                  'da': 'Damping parameters from 1D synthesis', \
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

            # contribution_1D
            elif label == 'bc':

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

            # Populations and departures
            elif label == 'bp' or label == 'bb':

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

                print(f'{label} found in 2 char')

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

                print(f'{label} found in 3 char')

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

                            # Valid class
                            return True

                    # 1.5D Contribution function
                    elif label == '2Dbc':

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

                            # Valid class
                            return True

                    # 1.5D height tau equal 1
                    elif label == '2Dbt':

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

                            # Valid class
                            return True

                    # Model atmosphere
                    elif label == '2Dat':

                        # Load 3D atmospheric model class
                        self.__object = _atmo_15D(self.__filename)

                        if self.__object is not None:

                            # Methods
                            self.get_filename = self.__get_filename
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

                            # Valid class
                            return True

                    # inversion_in
                    elif label == 'invi':

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

                            # Valid class
                            return True

                    # inversion_out
                    elif label == 'invo':

                        # Load inversion output class
                        self.__object = _inversion_out(self.__filename)

                        if self.__object is not None:

                            # Methods
                            self.get_filename = self.__get_filename
                            self.get_polarized = self.__get_polarized
                            self.get_nx = self.__get_nx
                            self.get_ny = self.__get_ny
                            self.get_nz = self.__get_nz
                            self.get_nl = self.__get_nl
                            self.get_nxy = self.__get_nxy
                            self.get_nxyz = self.__get_nxyz
                            self.get_dims = self.__get_dims
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

                            # Valid class
                            return True

                    # collisional rates and damping
                    elif label == '2Dct' or label == '2Dcl' or label == '2Dda':

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

                            # Valid class
                            return True

                    # Continuum quantities
                    elif label == '2Dba':

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

                            # Valid class
                            return True

                    # Populations and departures
                    elif label == '2Dbp' or label == '2Dbb':

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

                            # Valid class
                            return True

                    # Fail
                    else:

                        # Not valid class
                        return False

                        print(f'{label} found in 4 char')

                # Not any of the labels
                else:

                    print(f'{label} not found')
                    return False

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
                self.__verbose('You can call help(method=name,input=name) ' + \
                         'to get information about a method or input',True)

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
                                   f'{input}: {methods[method][0][input]}',True)

################################################################################

    # Parsers

    # stokes 1D, contribution 1D, tau 1D, jkq 1D, rkq 1D, stokes 15D,
    # contribution 15D, tau 15D, 3D atmos, inversion out, 15D cols,
    # 15D back, 15D popdep
    def __get_filename(self):
        return self.__object._get_filename()
    # 15D cols, 15D popdep
    def __get_type(self):
        return self.__object._get_type()
    # stokes 15D, inversion out
    def __get_polarized(self):
        return self.__object._get_polarized()
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
    # 3D atmos, inversion out, 15D back
    def __get_vars(self):
        return self.__object._get_vars()
    # 3D atmos, inversion out, 15D back
    def __get_vars_units(self):
        return self.__object._get_vars_units()
    # 3D atmos, 15D back
    def __get_vars_alias(self):
        return self.__object._get_vars_alias()
    # inversion out
    def __get_vars_atmo(self):
        return self.__object._get_vars_atmo()
    # inversion out
    def __get_vars_atmo_units(self):
        return self.__object._get_vars_atmo_units()
    # stokes 1D, contribution 1D, tau 1D, 1D popdep, stokes 15D,
    # contribution 15D, tau15D, inversion in, inversion out, 15D back
    def __get_nl(self):
        return self.__object._get_nl()
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
    # inversion out, 15D cols, 15D back, 15D popdep, 1D popdep
    def __get_dims(self):
        return self.__object._get_dims()
    # contribution 1D, jkq 1D, rkq 1D, 1D popdep, contribution 15D,
    # 3D atmos, inversion out, 15D cols, 15D back, 15D popdep
    def __get_nz(self):
        return self.__object._get_nz()
    # 15D cols, 15D popdep
    def __get_nentry(self):
        return self.__object._get_nentry()
    # jkq 1D, rkq 1D, 15D cols
    def __get_na(self):
        return self.__object._get_na()
    # jkq 1D
    def __get_nt(self):
        return self.__object._get_nt(atom=None)
    # inversion in
    def __get_los(self):
        return self.__object._get_los()
    # inversion in
    def __get_los_c(self,ix=None,iy=None):
        return self.__object._get_los_column(ix,iy)
    # inversion in
    def __get_los_p(self,ix=None,iy=None):
        return self.__object._get_los_plane()
    # stokes 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D
    def __get_th(self):
        return self.__object._get_th()
    # stokes 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D
    def __get_ph(self):
        return self.__object._get_ph()
    # stokes 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D
    def __get_mu(self):
        return self.__object._get_mu()
    # stokes 1D, contribution 1D, tau 1D, stokes 15D,
    # contribution 15D, tau15D, inversion in, inversion out, 15D back
    def __get_lambda(self,minl=None,maxl=None):
        return self.__object._get_lambda(minl,maxl)
    # contribution 1D, jkq 1D, rkq 1D
    def __get_height(self,minh=None,maxh=None):
        return self.__object._get_height(minh,maxh)
    # tau 1D
    def __get_height1d(self,minl=None,maxl=None):
        return self.__object._get_height(minl,maxl)
    # stokes 1D
    def __get_stokesi1d(self,minl=None,maxl=None):
        return self.__object._get_stokesi(minl,maxl)
    # stokes 1D
    def __get_stokesq1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesq(minl,maxl,fractional)
    # stokes 1D
    def __get_stokesu1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesu(minl,maxl,fractional)
    # stokes 1D
    def __get_stokesv1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesv(minl,maxl,fractional)
    # stokes 1D
    def __get_stokes1d(self,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokes(minl,maxl,fractional)
    # stokes 1D
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
        return self.__object._get_jkq(atom,transition,k,q,minh,maxh,sti)
    # rkq 1D
    def __get_rkq1d(self,atom=None,minh=None,maxh=None):
        return self.__object._get_rkq(atom,minh,maxh)
    # contribution 1D
    def __get_ctr1d(self,minl=None,maxl=None,minh=None,maxh=None):
        return self.__object._get_ctr(minl,maxl,minh,maxh)
    # 1D popdep
    def __get_column_1dcapd(self,ie=None):
        return self.__object._get_data(ie)
    # stokes 15D, inversion in
    def __get_stokesi15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_stokesi_column(ix,iy,minl,maxl)
    # stokes 15D, inversion in
    def __get_stokesq15d_c(self,ix=None,iy=None,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesq_column(ix,iy,minl,maxl,fractional)
    # stokes 15D, inversion in
    def __get_stokesu15d_c(self,ix=None,iy=None,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesu_column(ix,iy,minl,maxl,fractional)
    # stokes 15D, inversion in
    def __get_stokesv15d_c(self,ix=None,iy=None,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokesv_column(ix,iy,minl,maxl,fractional)
    # stokes 15D, inversion in
    def __get_stokes15d_c(self,ix=None,iy=None,minl=None,maxl=None,fractional=False):
        return self.__object._get_stokes_column(ix,iy,minl,maxl,fractional)
    # stokes 15D, inversion in
    def __get_linear15d_c(self,ix=None,iy=None,minl=None,maxl=None,fractional=False):
        return self.__object._get_linear_column(ix,iy,minl,maxl,fractional)
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
    def __get_stokesv15d_p(self,i=None,fractional=False):
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
    # tau 15D
    def __get_height15d_c(self,ix=None,iy=None,minl=None,maxl=None):
        return self.__object._get_column(ix,iy,minl,maxl)
    # tau 15D
    def __get_height15d_p(self,il=None):
        return self.__object._get_height_plane(il)
    # 3D atmo
    def __get_column_2dat(self,ix=None,iy=None,minh=None,maxh=None, \
                          mint=None,maxt=None,var=None):
        return self.__object._get_column(ix,iy,minh,maxh,mint,maxt,var)
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
        return self.__object._get_plane_node(var)
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

################################################################################
################################################################################
################################################################################
