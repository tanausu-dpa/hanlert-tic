from hanlertmodlib as hrtml
import numpy as np
#from astropy.io import fits
#from scipy.io import readsav
#from h5py

check = True

def main():
    ''' This routine helps building the binary data file for
        the 15D synthesis mode of HanleRT-TIC
    '''

    # Write the path to the output file
    output_file = ''

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
    helper = hrtml.atmo15D_class(output_file)

    ##################################################################
    ##################################################################
    # Here you need to provide the helper class with the data        #
    ##################################################################
    ##################################################################

    #
    # IMPORTANT! X is the slow axis, and Y is the fast axis. The
    # memory order is the important distinction, not the actual name
    #
    # IMPORTANT! If you only have a single spatial dimensions, you
    # still need to provide your arrays with two spatial dimensions,
    # so you must add a size 1 additional dimension for the other
    # spatial axis.
    #

    #
    # Precision
    ##################################################################

    #
    # The file can hold the model in single or double precision, so
    # its precision must be specified. If this is not specified here,
    # the precision will be assumed from the first variable that is
    # loaded.
    #
   #check = helper.set_precision('s')
    #
    # or
    #
   #check = helper.set_precision('d')
    #

    # Do not remove this sanity check
    if not check: sys.exit()

    #
    # Physical data
    ##################################################################

    #
    # We need to provide the 3D cubes as numpy arrays. If any or both
    # x or y are 1, the dimension must still exist with size 1.

    #
    # IMPORTANT: If the model is large, it is very recommended to
    # use numpy.memmap to handle the variables

    #
    # We can inject the variables one by one
    #
    # Things to take into consideration:
    #
    #  - It is necessary to specify either heights or optical depths.
    #    IMPORTANT! It is tau, not log_10(tau) !!!
    #
    #  - It is necessary to specify either gas pressure (Pg),
    #    electron pressure (Pe-), electron density (ne-), or
    #    electron density (ne-) and hydrogen number densities (nH_*).
    #    More than one can be specified, the one used is configured
    #    in the input file of HanleRT-TIC.
    #
    #  - Only relevant if writing the hydrogen minus number density.
    #    The model distinguishes between total hydrogen number
    #    density, atomic hydrogen number density, and hydrogen
    #    minus number density. The atomic number density includes
    #    only the bound (neutral) and free (positive) states. The
    #    total includes the bound (negative) and molecular. If
    #    you do not have explicit values for atomic or total, they
    #    will be automatically calculated from the other inputs.
    #
    # Now we load the variables. It is typical that heights (or
    # sometimes optical depth) are the same for all pixels. If
    # this is the case, it is easier to broadcast them. Here we
    # assume nx, ny, and nz are the sizes already specified.
    #
   #check = helper.load_data(np.broadcast_to(my_h,(nx,ny,nz)), \
   #                         'h',ix=0,iy=1)
    # or
   #check = helper.load_data(np.broadcast_to(my_tau,(nx,ny,nz)), \
   #                         'tau',ix=0,iy=1)
    #
    # If we have them in 3D, then we just plug them as any other
    # variable (see below). In the following commands, the three
    # index positions are specified, but only two are necessary.
    # The list of available variables is:
    #
    #
    # Height. To be specified in km
   #check = helper.load_data(my_h,'h',ix=<index_x_dim>, \
   #                                  iy=<index_y_dim>, \
   #                                  iz=<index_z_dim>)
   #
    # Optical depth. It is dimensionless
   #check = helper.load_data(my_tau,'tau',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Temperature. To be specified in kelvin (K)
   #check = helper.load_data(my_T,'T',ix=<index_x_dim>, \
   #                                  iy=<index_y_dim>, \
   #                                  iz=<index_z_dim>)
   #
    # Gas pressure. To be specified in cgs (dyn cm^-2)
   #check = helper.load_data(my_Pg,'Pg',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Density. To be specified in cgs (g cm^-3)
   #check = helper.load_data(my_rho,'rho',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Magnetic field X component. To be specified in gauss (G)
   #check = helper.load_data(my_Bx,'Bx',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Magnetic field Y component. To be specified in gauss (G)
   #check = helper.load_data(my_By,'By',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Magnetic field Z component. To be specified in gauss (G)
   #check = helper.load_data(my_Bz,'Bz',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Velocity X component. To be specified in km s^-1
   #check = helper.load_data(my_vx,'vx',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Velocity Y component. To be specified in km s^-1
   #check = helper.load_data(my_vy,'vy',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Velocity Z component. To be specified in km s^-1
   #check = helper.load_data(my_vz,'vz',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Microturbulent velocity. To be specified in km s^-1
   #check = helper.load_data(my_vm,'vmi',ix=<index_x_dim>, \
   #                                     iy=<index_y_dim>, \
   #                                     iz=<index_z_dim>)
   #
    # Electron pressure. To be specified in cgs (dyn cm^-2)
   #check = helper.load_data(my_Pe,'Pe',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Electron number density. To be specified in cgs (cm^-3)
   #check = helper.load_data(my_ne,'ne',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)
   #
    # Total hydrogen number density. To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nht,'nht',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Atomic hydrogen number density. To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nha,'nha',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Hydrogen minus number density. To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nhm,'nhm',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Hydrogen ground level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh0,'nh0',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Hydrogen first excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh1,'nh1',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Hydrogen second excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh2,'nh2',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Hydrogen third excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh3,'nh3',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Hydrogen fourth excited level number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_nh4,'nh4',ix=<index_x_dim>, \
   #                                      iy=<index_y_dim>, \
   #                                      iz=<index_z_dim>)
   #
    # Proton number density.
    # To be specified in cgs (cm^-3)
   #check = helper.load_data(my_np,'np',ix=<index_x_dim>, \
   #                                    iy=<index_y_dim>, \
   #                                    iz=<index_z_dim>)


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
