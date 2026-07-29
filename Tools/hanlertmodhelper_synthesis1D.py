import hanlertmodlib as hrtml
import numpy as np
#from astropy.io import fits
#from scipy.io import readsav
#from h5py

check = True

def main():
    ''' This routine helps building the binary data file for
        the 1D synthesis mode of HanleRT-TIC
    '''

    # Write the path to the output file for the model atmosphere
    atmosphere_file = ''
    # and for the magnetic field
    magnetic_file = ''

    ##################################################################
    ##################################################################
    # Here you should read your data declaring some unique variables #
    ##################################################################
    ##################################################################

    # To read fits files, you need the astropy library. Remove the
    # comment in the header of this file for the fits import.
    # Instructions to read fits with astropy can be found at
    # https://docs.astropy.org/en/stable/io/fits/


    # To read IDL save files, you need the scipy library. Remove
    # the comment in the header of this file for the sav import.
    # To read:
    # variable = readsav(<path_to_file>)
    # This creates a dictionary, you can check the content with
    # print(list(variable))


    # To read hdf5 files you need the h5py libreary. Remove the
    # comment in the header of this file for the hdf5 import.
    # The documentation for h5py is hosted at
    # https://docs.h5py.org/en/stable/index.html


    ##################################################################
    ##################################################################

    # Initialize the instance
    # The magnetic field file is optional
    helper = hrtml.atmo1D_class(atmosphere_file,magnetic_file)

    ##################################################################
    ##################################################################
    # Here you need to provide the helper class with the data        #
    ##################################################################
    ##################################################################

    #
    # Comment 
    ##################################################################

    #
    # The 1D file is in ASCII and therefore can include a header
    # with a comment. 
    #
   #check = helper.set_comment('My insightful comments')

    #
    # We can also specify a name for the model
    #
   #check = helper.set_name('COOL MODEL')

    # Do not remove this sanity check
    if not check: sys.exit()


    #
    # Wavelength
    ##################################################################

    #
    # If an optical depth is specified (or is desired as output of
    # HanleRT-TIC), the continuum wavelength must be specified. By
    # default it is 500 nm.
    #
   #check = helper.set_wavelength(wavelength_in_nm)

    # Do not remove this sanity check
    if not check: sys.exit()

    #
    # Physical data
    ##################################################################

    #
    # We need to provide the 1D numpy arrays.

    #
    # We can inject the variables one by one
    #
    # Things to take into consideration:
    #
    #  - Only height or optical depth can be specified, not both.
    #    IMPORTANT! It is tau, not log_10(tau) !!!
    #
    #  - Only one of gas pressure, electron pressure, density, or
    #    electron number density can be specified. If electron
    #    number density is specified, then hydrogen number
    #    densities can be specified as well.
    #
    # Now we load the variables
    #
   #check = helper.load_data(my_h,'h')
    # or
   #check = helper.load_data(my_tau,'tau')
    #
    # Temperature. To be specified in kelvin (K)
   #check = helper.load_data(my_T,'T')
   #
    # Gas pressure. To be specified in cgs (dyn cm^-2)
   #check = helper.load_data(my_Pg,'Pg')
    # or density. To be specified in cgs (g cm^-3)
   #check = helper.load_data(my_rho,'rho')
    # or electron pressure. To be specified in cgs (dyn cm^-2)
   #check = helper.load_data(my_Pe,'Pe')
    # or electron density. To be specified in cgs (g cm^-3)
   #check = helper.load_data(my_rhoe,'rhoe')
    # or electron number density. To be specified in cgs (cm^-3)
   #check = helper.load_data(my_ne,'ne')
    #
    # Magnetic field X component. To be specified in gauss (G)
   #check = helper.load_data(my_Bx,'Bx')
   #
    # Magnetic field Y component. To be specified in gauss (G)
   #check = helper.load_data(my_By,'By')
   #
    # Magnetic field Z component. To be specified in gauss (G)
   #check = helper.load_data(my_Bz,'Bz')
   #
    # Velocity X component. To be specified in km s^-1
   #check = helper.load_data(my_vx,'vx')
   #
    # Velocity Y component. To be specified in km s^-1
   #check = helper.load_data(my_vy,'vy')
   #
    # Velocity Z component. To be specified in km s^-1
   #check = helper.load_data(my_vz,'vz')
   #
    # Microturbulent velocity. To be specified in km s^-1
   #check = helper.load_data(my_vm,'vmi')
   #
    # Hydrogen ground level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh0,'nh0')
   #
    # Hydrogen first excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh1,'nh1')
   #
    # Hydrogen second excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh2,'nh2')
   #
    # Hydrogen third excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh3,'nh3')
   #
    # Hydrogen fourth excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh4,'nh4')
   #
    # Proton number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_np,'np')

    # Do not remove this sanity check
    if not check: sys.exit()

    #
    # We try to create the file here
    ##################################################################

    check = helper.create_file()

    # Do not remove this sanity check
    if not check: sys.exit()

    ##################################################################
    ##################################################################

    print('Your model atmosphere file is ' + \
          'ready at {0}'.format(output_file))


if __name__ == '__main__':
    main()
