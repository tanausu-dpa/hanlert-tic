The magnetic field files have the following structure line by line.

1. Number of height nodes. This number must be the same that the
   number of nodes in the :ref:`model atmosphere <atmo_file>`.

   ::

     70


   .. note::

      If the specified number of nodes is exactly 1, then the
      specified magnetic field will be replicated in all heights,
      that is, to specify an homogeneous magnetic field it is
      enough to indicate a single height node in the magnetic field
      file.

   .. tip::

      Note that in order to define an homogeneous magnetic field
      it is not necessary to specify a magnetic field file, see
      **BFIELD_INPUT** in :ref:`control file <inputs_control>` for
      further details.

2. The units of the angles in the file. Admitted strings are,

   * *RAD*. Angles will be in radians.
   * *DEG*. Angles will be in degrees.

   ::

      DEG

3. As many lines as number of height nodes with three floating
   point numbers indicating:

   #. The magnetic field strength in gauss.
   #. The angle of the magnetic field vector with respect to the vertical
      (in radians if *RAD* and in degrees if *DEG*).
   #. The angle of the magnetic field vector projection on the plane
      normal to the vertical with respect to the reference drection of
      zero azimuth (in radians if *RAD* and in degrees if *DEG*).

   ::

      10.0 15.0 0.0
      12.3 16.4 0.0
        .    .   .
        .    .   .
        .    .   .
