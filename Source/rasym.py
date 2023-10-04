import sys, math, os

#####################
# razym()
#
# Tanaus\'u del Pino Alem\'an
#
# 06/29/2022(US):  V3.0.0 - Changed global version (TdPA)
#
# 03/17/2021(US):  V2.0.0 - Changed global version (TdPA)
#
# 12/16/2020(US):  V1.0.0 - First version
#
#####################


def rasym():
  ''' Reads the magnetic field file specified as argument.
  '''

  def abort(f,name):
    f.close()
    f = open(name,'w')
    f.write('-1')
    f.close()
    sys.exit()

  def verbose(msg, fil, verb):

    # If being verbose
    if (verb):
      print((msg+' in rasym.py'))
    else:
      exist = os.path.isfile(fil)
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      fv.write(msg+' in rasym.py\n')
      fv.close()

  # Argument control
  if len(sys.argv) < 1:
   #sys.exit(' # At least one argument needed')
    sys.exit()
  try:
    dni = sys.argv[2]
  except:
    dni = '000000000'
  if len(sys.argv) > 3:
    verbosity = False
    verbfile = sys.argv[3]
  else:
    verbosity = True
    verbfile = ''

  try:
    f=open(sys.argv[1],'r')
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

  # Output file
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
    if iline == len(lines):
        break

    # Components
    val = lines[iline].split()

    # If four elements
    if len(val) == 4:

      try:

        K = int(val[0])
        Q = int(val[1])
        nz = 1
        ar = float(val[2])
        ai = float(val[3])

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

      except:
        verbose(' # Error reading K, Q, and factor from ' + \
                '{0},{1},{2},{3}'.format(*val), \
                verbfile,verbosity)
        abort(f, filename)

    # If three elements
    elif len(val) == 3:

      try:

        K = int(val[0])
        Q = int(val[1])
        nz = int(val[2])

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

      except:
        verbose(' # Error reading K, Q, and nz from ' + \
                '{0},{1},{2}'.format(*val), \
                verbfile,verbosity)
        abort(f, filename)

      try:

        # For each nz
        for iz in range(nz):

          # Advance line
          iline += 1

          # If EoF
          if iline == len(lines):
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

          if (ar*ar + ai*ai) < 0. or math.sqrt(ar*ar + ai*ai) > 1.:
            verbose(' # Input asymmetry absolute value must ' + \
                    'be in the (0,1] range, your input is ' + \
                    '{0} + i{1}'.format(ar,ai), verbfile, verbosity)
            abort(f, filename)

          # Add to writing buffer
          out.append('{0} {1}'.format(ar,ai))

      except:
        verbose(' # Error reading anisotropy factor from line ' + \
                '{0} {1}'.format(*lval), \
                verbfile,verbosity)
        abort(f, filename)

  # Finish writing file
  f.write('{0}\n'.format(nentry))
  for o in out:
      f.write(o+'\n')
  f.close()

if __name__ == "__main__":
    rasym()
