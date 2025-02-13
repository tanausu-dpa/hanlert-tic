import sys, math, os

#####################
# rBField()
#
# Tanaus\'u del Pino Alem\'an (IAC)
#
# 17/12/2024:  V4.0.0 - Changed global version (TdPA)
#
#####################


def rBField():
  ''' Reads the magnetic field file specified as argument.
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
  def verbose(msg, fil, verb):
    # If being verbose
    if (verb):
      # Just print
      print((msg+' in rbfield.py'))
    else:
      # Check file exists
      exist = os.path.isfile(fil)
      # Open to write or append
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      # Write in file and close
      fv.write(msg+' in rbfield.py\n')
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
    verbose(' # No field file found', verbfile, verbosity)
    filename = 'tmp_bfield_'+dni
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
  filename='tmp_bfield_'+dni
  f=open(filename,'w')
  f.write('1\n')

  #
  # File build
  #

  # Get dimenions
  val = lines[0].split()[0]
  try:
    NZ=int(float(val))
  except:
    verbose(' # First line must be a number', \
            verbfile, verbosity)
    abort(f, filename)
  f.write('{0:7d}\n'.format(NZ))
  if(NZ != (len(lines)-2)):
    verbose(' # Bad dimensions, number of lines in ' + \
            'magnetic file is not equal to number ' + \
            'specified', verbfile, verbosity)
    abort(f, filename)

  # Get angular units
  val = lines[1].split()[0]
  try:
    val = val.upper()
  except:
    verbose(' # Second line in magnetic file must be a string', \
            verbfile, verbosity)
    abort(f, filename)
  if val == 'RAD':
    f.write('RAD\n')
  elif val == 'DEG':
    f.write('DEG\n')
  else:
    verbose(' # Second line in magnetic file not understood,' + \
            ' set to DEG(DEFAULT)', verbfile, verbosity)
    f.write('DEG\n')

  # For each height
  for i in range(NZ):
    # Read field in polar coordinates
    val = lines[2+i].split()
    if(len(val) != 3):
      verbose(' # The format must have three columns', \
              verbfile, verbosity)
      abort(f, filename)
    for j in range(3):
      f.write('{0:22.16e} '.format(float(val[j])))
    f.write('\n')

  # Close file
  f.close()

if __name__ == "__main__":
    rBField()
