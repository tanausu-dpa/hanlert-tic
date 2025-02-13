import sys, math, os

#####################
# rasym()
#
# Tanaus\'u del Pino Alem\'an (IAC)
#
# 13/12/2024:  V4.0.0 - Changed global version (TdPA)
#
#####################


def rasym():
  ''' Reads the ad-hoc radiation field tensors file specified as
      argument
  '''

  # Aborting method
  def abort(f,name):
    # Close file
    f.close()
    # Reset file and just write -1 to flag failure
    f = open(name,'w')
    f.write('-1')
    f.close()
    # Leave
    sys.exit()

  # Verbose routine
  def verbose(msg,fil,verb):
    # If being verbose
    if (verb):
      # Just print
      print((msg+' in rasym.py'))
    else:
      # Check file exists
      exist = os.path.isfile(fil)
      # Open to write or append
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      # Write in file and close
      fv.write(msg+' in rasym.py\n')
      fv.close()

  #
  # Argument control
  #

  # Requires one argument
  if len(sys.argv) < 1:
    sys.exit()

  # Try getting ID
  try:
    dni = sys.argv[2]
  except:
    dni = '000000000'

  # If more arguments, there is verbosity file
  if len(sys.argv) > 3:
    verbosity = False
    verbfile = sys.argv[3]
  else:
    verbosity = True
    verbfile = ''

  # Try to open file
  try:
    f=open(sys.argv[1],'r')
  # Failed to open file
  except:
    verbose(' # No asymmetry file found', verbfile, verbosity)
    filename = 'tmp_asym_'+dni
    f = open(filename,'w')
    abort(f, filename)

  # Read file
  lines=list(f)
  f.close()

  # Manage Comments
  lines_n = []
  for line in lines:
    line = line.strip()
    j = line.find('!')
    if j != -1:
      line = line[:j]
    if len(line.split()) == 0:
      continue
    else:
      lines_n.append(line.upper())
  lines = lines_n

  # Start output file
  filename='tmp_asym_'+dni
  f=open(filename,'w')
  f.write('1\n')

  # Initialize entries counter, status, and buffer for writing
  nentry = 0
  ninentry = 0
  out = []
  iline = -1

  # For each line
  while True:

    # Advance line
    iline += 1

    # EoF
    if iline == len(lines):
        break

    # Components
    val = lines[iline].split()

    # If four elements
    if len(val) == 4:

      # Try parsing
      try:

        # Parse
        K = int(val[0])
        Q = int(val[1])
        nz = 1
        ar = float(val[2])
        ai = float(val[3])

        # Sanity check
        if K < 1 or K > 2:
          verbose(' # Input asymmetry must have ' + \
                  'multipole K in the [1,2] range, '+ \
                  'your input is {0}'.format(K), \
                  verbfile,verbosity)
          abort(f, filename)
        if Q < 0 or Q > K:
          verbose(' # Input asymmetry must have ' + \
                  'multipole Q in the [0,K] range, ' + \
                  'your input is {0}'.format(Q), \
                  verbfile,verbosity)
          abort(f, filename)
        if (ar*ar + ai*ai) <= 0. or math.sqrt(ar*ar + ai*ai) > 1.:
          verbose(' # Input asymmetry absolute value must ' + \
                  'be in the (0,1] range, your input is ' + \
                  '{0} + i{1}'.format(ar,ai), verbfile, verbosity)
          abort(f, filename)

        # Add entry
        nentry += 1
        out.append('{0} {1} {2}'.format(K,Q,nz))
        out.append('{0} {1}'.format(ar,ai))

      # Could not parse
      except:
        # Interpretation error
        verbose(' # Error reading K, Q, and factor from ' + \
                '{0},{1},{2},{3}'.format(*val), \
                verbfile,verbosity)
        abort(f, filename)

    # If three elements
    elif len(val) == 3:

      # Try parsing
      try:

        # Parse
        K = int(val[0])
        Q = int(val[1])
        nz = int(val[2])

        # Sanity check
        if K < 1 or K > 2:
          verbose(' # Input asymmetry must have ' + \
                  'multipole K in the [1,2] range, '+ \
                  'your input is {0}'.format(K), \
                  verbfile,verbosity)
          abort(f, filename)
        if Q < 0 or Q > K:
          verbose(' # Input asymmetry must have ' + \
                  'multipole Q in the [0,K] range, ' + \
                  'your input is {0}'.format(Q), \
                  verbfile,verbosity)
          abort(f, filename)
        if nz < 3 and nz != 1:
          verbose(' # Input node number must be larger ' + \
                  'than 3, your is {0}'.format(nz), \
                  verbfile,verbosity)
          abort(f, filename)

        # Add entry
        nentry += 1
        out.append('{0} {1} {2}'.format(K,Q,nz))

      # Could not parse
      except:
        # Interpretation error
        verbose(' # Error reading K, Q, and nz from ' + \
                '{0},{1},{2}'.format(*val), \
                verbfile,verbosity)
        abort(f, filename)

      # Try getting stratification
      try:

        # For each nz
        for iz in range(nz):

          # Advance line
          iline += 1

          # If EoF
          if iline == len(lines):
            # Error
            verbose(' # End of file {0}'.format(filename) + \
                    'reached before finishing entry', \
                    verbfile,verbosity)
            abort(f, filename)

          # Get components
          lval = lines[iline].split()

          # Check columns
          if len(lval) != 2:
            verbose(' # The following line has more than two ' + \
                    'values and should be an entry' + lines[iline], \
                    verbfile,verbosity)
            abort(f, filename)

          # Factors
          ar = float(lval[0])
          ai = float(lval[1])

          # Sanity check
          if (ar*ar + ai*ai) < 0. or math.sqrt(ar*ar + ai*ai) > 1.:
            verbose(' # Input asymmetry absolute value must ' + \
                    'be in the (0,1] range, your input is ' + \
                    '{0} + i{1}'.format(ar,ai), verbfile, verbosity)
            abort(f, filename)

          # Add to writing buffer
          out.append('{0} {1}'.format(ar,ai))

      # Could not get stratification
      except:
        # Error
        verbose(' # Error reading anisotropy factor from line ' + \
                '{0} {1}'.format(*lval), \
                verbfile,verbosity)
        abort(f, filename)

  # Finish writing file
  f.write('{0}\n'.format(nentry))
  for o in out:
      f.write(o+'\n')

  # And close
  f.close()

if __name__ == "__main__":
    rasym()
