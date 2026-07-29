import hanlertmodlib as hrtml
import numpy as np
#from astropy.io import fits
#from scipy.io import readsav
#from h5py

check = True

def main():
    ''' This routine helps building the binary data file for
        the inversion mode of HanleRT-TIC
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
    helper = hrtml.invi_class(output_file)

    ##################################################################
    ##################################################################
    # Here you need to provide the helper class with the data        #
    ##################################################################
    ##################################################################

    #
    # IMPORTANT! X is the slow axis, and Y is the fast axis. The
    # memory order is the important distinction, not the actual name.
    #
    # IMPORTANT! If you only have a single spatial dimensions, you
    # still need to provide your arrays with two spatial dimensions,
    # so you must add a size 1 additional dimension for the other
    # spatial axis.
    #

    #
    # Wavelength axis
    ##################################################################

    #
    # Here you would load or create a numpy array, dubbed
    # 'my_wavelength' in the example, in nanometers.
   #my_wavelength = TODO define wavelengths

    # If you wavelengths are in air, use the following to convert to
    # vacuum (assuming wavelengths in nanometers).
   #my_wavelength = hrtml.airtovacuum(my_wavelength*1e1)*1e-1

    #
    # To load the wavelengths. They must be in nanometers.
   #check = helper.load_wavelength(my_wavelength)

    # Do not remove this sanity check
    if not check: sys.exit()

    #
    # Units
    ##################################################################

    #
    # The code expects everything to be in the international system
    # of units, in units of hertz (Hz). This applies to the intensity,
    # the polarization Stokes parameters, the sigma (errors), and
    # the diffuse light.
    #
    # If your variables are in wavelength units, intead of frequency
    # units, the transformation corresponds to multiply by:
    #
    #    (my_wavelength_u)*(my_wavelenth*1e-7)/2.99792458e10
    #
    # where my_wavelength is the wavelength axis in nanometers that
    # was loaded before, and my_wavelength_u is the same axis, with
    # the same units that the wavelength band of your units. The
    # wavelength band is typically in Angstrom, so it would be
    #     my_wavelength_u = my_wavelength*1e1
    #
    # If your variables are in cgs units, in units of frequency, to
    # transform into SI it is enough to multiply by 1e-3
    #
    # If your variables are in units of the continuum, you can
    # use the included hanlertmodallen.py to get the continuum
    # intensity tabulated by Allen (or any other method you see fit
    # to calibrate in intensity). To get the intensity at a certain
    # wavelength:
    #    import hanlertmodallen as hrta
    #    allen = hrta.allen_class()
    #    Ic = allen.get_radiation(wavelength_in_nm)*1e-3
    #


    #
    # Stokes data
    ##################################################################

    #
    # Here we would define or load the Stokes parameters. In the
    # example they are stored in the 'my_data' variable

    #
    # To load your data cube. They must be in international system of
    # units. For different dimensions of the array:
    #
    #    - 1D: It is assumed that it is a single intensity spectrum,
    #          so the dimension is automatically wavelength.
    #
   #check = helper.load_data(my_data)
    #
    #    - 2D: It is assumed that it is a single polarized spectrum.
    #          The helper tries to find the Stokes dimension, but you
    #          can provide it.
    #
   #check = helper.load_data(my_data)
    # or
   #check = helper.load_data(my_data,istk=<index_stokes_dim>)
    # or
   #check = helper.load_data(my_data,istk=<index_stokes_dim>, \
   #                                 il=<index_lambda_dim>)
    #
    #    - 3D: It is assumed that it is an only intensity data cube.
    #          You need to provide at least two of the three
    #          indexes indicating the position of the dimension in
    #          your array.
    #
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 iy=<index_y_dim>)
    # or
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 il=<index_lambda_dim>)
    # or
   #check = helper.load_data(my_data,iy=<index_y_dim>, \
   #                                 il=<index_lambda_dim>)
    # or
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 iy=<index_y_dim>, \
   #                                 il=<index_lambda_dim>)
    #
    #    - 4D: It is assumed that it is a full Stokes data cube.
    #          You need to provide at least three of the four
    #          indexes indicating the position of the dimension in
    #          your array.
    #
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 iy=<index_y_dim>, \
   #                                 il=<index_lambda_dim>)
    # or
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 iy=<index_y_dim>, \
   #                                 istk=<index_stokes_dim>)
    # or
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 il=<index_lambda_dim>, \
   #                                 istk=<index_stokes_dim>)
    # or
   #check = helper.load_data(my_data,iy=<index_y_dim>, \
   #                                 il=<index_lambda_dim>, \
   #                                 istk=<index_stokes_dim>)
    # or
   #check = helper.load_data(my_data,ix=<index_x_dim>, \
   #                                 iy=<index_y_dim>, \
   #                                 il=<index_lambda_dim>, \
   #                                 istk=<index_stokes_dim>)
    #

    # Do not remove this sanity check
    if not check: sys.exit()

    #
    # Line of sight data
    ##################################################################

    #
    # The line of sight must be specified as the heliocentric and
    # azimuthal angles (only really critical for scattering
    # polarization, as this defines the reference direction for
    # positive Stokes Q) in radians
    #
    # If the LOS is shared for the whole field of view, a one
    # dimensions size 2 array is to be provided, e.g.,
    # my_los = np.array([heliocentric,azimuth])
    #
   #check = helper.load_los(my_los)
    #
    # If the LOS changes pixel by pixel, the array must have three
    # dimensions, and the indexes for the x and y dimensions need to
    # be specified.
    #
   #check = helper.load_los(my_los,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>)
    #

    # Do not remove this sanity check
    if not check: sys.exit()


    #
    # Sigma data
    ##################################################################

    #
    # To load the sigma for your data (errors). They must be in
    # international system of units. For different types of data:
    #
    # o Only intensity data:
    #
    #    - A single number assumes that the same error is to be
    #      considered for every wavelength and pixel
    #
   #check = helper.load_sig(my_sig,polarization=False)
    #
    #    - A 1D array assumes that the same error is to be considered
    #      for every pixel, but it is wavelength dependent
    #
   #check = helper.load_sig(my_sig,polarization=False)
    #
    #    - A 2D array assumes that the error is not wavelength
    #      dependent, but changes for each pixel. You must specify the
    #      index of at least 1 spatial dimension
    #
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               polarization=False)
    # or
   #check = helper.load_sig(my_sig,iy=<index_y_dim>, \
   #                               polarization=False)
    # or
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               polarization=False)
    #
    #    - A 3D array assumes that the error is both wavelength and
    #      pixel dependent. You must specify at least two indexes
    #      for the spatial and wavelength dimensions.
    #
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               polarization=False)
    # or
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               polarization=False)
    # or
   #check = helper.load_sig(my_sig,iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               polarization=False)
    #
    # o Polarization data:
    #
    #    - A 1D array assumes that the same error is to be
    #      considered for every wavelength and pixel. This
    #      dimension must be size 4 (Stokes parameters)
    #
   #check = helper.load_sig(my_sig,polarization=True)
    #
    #    - A 2D array assumes that the same error is to be
    #      considered for every pixel, but it is wavelength
    #      dependent. The helper tries to find the Stokes dimension,
    #      but you can provide it.
    #
   #check = helper.load_sig(my_sig,polarization=True)
    # or
   #check = helper.load_sig(my_sig,il=<index_lambda_dim>, \
   #                               polarization=True)
    #
    #    - A 3D array assumes that the same error is not wavelength
    #      dependent, but changes for each pixel. You must specify at
    #      least two indexes for the spatial and Stokes dimensions.
    #
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_sig(my_sig,iy=<index_y_dim>, \
   #                               istk=<index_Stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               istk=<index_Stokes_dim>, \
   #                               polarization=True)
    #
    #    - A 4D array assumes that the error is both wavelength and
    #      pixel dependent. You must specify at least three indexes
    #      for the spatial, wavelength, and Stokes dimensions.
    #
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_sig(my_sig,iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_sig(my_sig,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)


    # Do not remove this sanity check
    if not check: sys.exit()


    #
    # Diffuse light data
    ##################################################################

    #
    # To load the diffuse light for your data. They must be in
    # international system of units. For different types of data:
    #
    # o Only intensity diffuse light:
    #
    #    - A 1D array assumes that the same profile for every pixel.
    #
   #check = helper.load_dif(my_dif,polarization=False)
    #
    #    - A 2D array assumes that the profile changes pixel by
    #      pixel. You must specify the index of at least 1 spatial
    #      dimension
    #
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               polarization=False)
    # or
   #check = helper.load_dif(my_dif,iy=<index_y_dim>, \
   #                               polarization=False)
    # or
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               polarization=False)
    #
    # o Polarized diffuse light:
    #
    #    - A 3D array assumes the same profile for every pixel.
    #      You must specify at least two indexes for the spatial
    #      and Stokes dimensions.
    #
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_dif(my_dif,iy=<index_y_dim>, \
   #                               istk=<index_Stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               istk=<index_Stokes_dim>, \
   #                               polarization=True)
    #
    #    - A 4D array assumes that the profile changes pixel by
    #      pixel. You must specify at least three indexes for the
    #      spatial, wavelength, and Stokes dimensions.
    #
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_dif(my_dif,iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)
    # or
   #check = helper.load_dif(my_dif,ix=<index_x_dim>, \
   #                               iy=<index_y_dim>, \
   #                               il=<index_lambda_dim>, \
   #                               istk=<index_stokes_dim>, \
   #                               polarization=True)


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

    print('Your data file is ready at {0}'.format(output_file))


if __name__ == '__main__':
    main()
