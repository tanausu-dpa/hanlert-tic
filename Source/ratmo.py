import sys, math, os

#####################
# rAtmo()
#
# Tanaus\'u del Pino Alem\'an (IAC)
#
# 17/12/2024:  V4.0.0 - Changed global version (TdPA)
#
#####################

def rAtmo():
  ''' Reads the input atmospheric file specified as argument.
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
      print((msg+' in ratmo.py'))
    else:
      # Check file exists
      exist = os.path.isfile(fil)
      # Open to write or append
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      # Write in file and close
      fv.write(msg+' in ratmo.py\n')
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
    verbose(' # No atmospheric file found', verbfile, verbosity)
    filename = 'tmp_atmo_'+dni
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
    j = line.find('*')
    if j != -1:
      line = line[:j]
    if len(line.split()) == 0:
      continue
    else:
      lines_n.append(line)
  lines = lines_n

  # Start output file
  filename = 'tmp_atmo_'+dni
  f = open(filename,'w')
  f.write('1\n')

  #
  # File build
  #

  # Get model's name
  iline = 0
  line = lines[iline].strip()
  name = line

  # Type of scale
  iline += 1
  line = lines[iline].strip()
  scale = line.split()
  # If info on reference frequency
  if len(scale) > 2:
    try:
      tf = '{0}\n'.format(1e2/float(scale[2]))
    except ValueError:
      verbose(' # Scale wavelength seems to not be float', \
              verbfile, verbosity)
      abort(f,filename)
    except:
      msg = ' # Critical unknown problem with reference ' + \
              'wavelength\n'
      for err in sys.exc_info()[:2]:
          msg += '   ' + err + '\n'
      verbose(msg,verbfile, verbosity)
      abort(f,filename)
  # If no info on reference frequency, assume 500 nm (0.2 kaiser)
  else:
    tf = '0.2\n'
  # Write scale
  if 'HEIGHT' in scale[0]:
    f.write('H\n')
  elif 'TAU' in scale[0]:
    f.write('T\n')
  else:
    verbose(' # Scale not recognized in model atmosphere', \
            verbfile, verbosity)
    abort(f,filename)
  # Write reference frequency
  f.write(tf)

  # Log(g), only used in inversion mode
  iline += 1
  line = lines[iline].strip()
  try:
    lst = list(line.lower())
    while lst.count('d') > 0:
        lst[lst.index('d')] = 'e'
    line = ''.join(lst)
    logg = float(line)
  except:
    verbose(' # Gravity in atmosphere must be a number', \
            verbfle, verbosity)
    abort(f,filename)
  f.write(line+'\n')

  # Number of height nodes
  iline += 1
  line = lines[iline]
  try:
    lst = list(line.lower())
    while lst.count('d') > 0:
        lst[lst.index('d')] = 'e'
    line = ''.join(lst)
    NZ = int(line.strip())
  except:
    verbose(' # Nodes in atmosphere must be number', \
            verbfile, verbosity)
    abort(f,filename)
  f.write(line+'\n')

  # Initialize atmo to write
  thermo = []
  dens = []
  nd = False

  #
  # Thermodynamical quantities
  #

  # For each height
  for iz in range(NZ):
    iline += 1
    line = lines[iline].strip()
    columns = line.split()
    # 5 numbers are required
    if len(columns) < 5:
      verbose(' # Needed at least five columns in ' + \
              'data block 1, got {0}'.format(len(colums)), \
              verbfile, verbosity)
      abort(f,filename)
    # If less than 7, fill horizontal velocity components with 0
    while len(columns) < 7:
      columns.append('0')
      line += ' 0'
    # Check they are numbers
    for ii in range(len(columns)):
      try:
        lst = list(columns[ii].lower())
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        columns[ii] = ''.join(lst)
        columns[ii] = float(columns[ii])
      except:
        verbose(' # Non numerical value in atmosphere', \
                verbfile, verbosity)
        abort(f,filename)
    # Add to buffer to write
    thermo.append(line+'\n')

  # Read one more line and check if it is an specifier
  iline += 1
  line = lines[iline].strip().lower()

  # If specifying only electron number density
  if 'ne' in line:
    f.write('1\n')

  # If specifying only electron pressure
  elif 'pe' in line:
    f.write('2\n')

  # If specifying only electron mass density
  elif 'rhoe' in line:
    f.write('3\n')

  # If specifying only gas pressure
  elif 'pg' in line:
    f.write('4\n')

  # If specifying only gas density
  elif 'rho' in line:
    f.write('5\n')

  # If we expect number densities
  else:

    # Flag number densities
    f.write('0\n')
    nd = True

    # Go back one line
    iline -= 1

    #
    # Hydrogen densities
    #

    # For each height
    for iz in range(NZ):
      iline += 1
      line = lines[iline].strip()
      columns = line.split()
      # Required at least two numbers
      if (len(columns)) < 2:
        verbose(' # Needed at least two columns in ' + \
                'data block 2, got {0}'.format(len(columns)), \
                verbfile, verbosity)
        abort(f,filename)
      changed = False
      # If less than 6 numbers, fill with zero the populations
      # of the missing levels until 5 neutral levels
      while (len(columns)) < 6:
        columns.insert(-1,'0')
        changed = True
      # Reconstruct the line if we changed it
      if changed:
        line = ' '.join(columns)
      # Check the line is made of numbers
      for ii in range(len(columns)):
        try:
          lst = list(columns[ii].lower())
          while lst.count('d') > 0:
            lst[lst.index('d')] = 'e'
          columns[ii] = ''.join(lst)
          columns[ii] = float(columns[ii])
        except:
          verbose(' # Non numerical value in atmosphere', \
                  verbfile, verbosity)
          abort(f,filename)
      dens.append(line+'\n')

  # Output thermodynamic quantities
  for line in thermo:
    f.write(line)

  # Output number densities
  for line in dens:
    f.write(line)

  # If number densities
  if nd:

    # Try to read helium number density
    try:
      outhe = ['1']
      # For each height
      for iz in range(NZ):
        iline += 1
        line = lines[iline].strip()
        columns = line.split()
        # Require at least three columns
        if (len(columns)) < 3:
          raise Exception
        changed = False
        # If less than for, add a zero in neutral stage
        while (len(columns)) < 4:
          columns.insert(1,'0')
          changed = True
        # If we added zeros, reconstruct the line
        if changed:
          line = ' '.join(columns)
        # Check the line is made of numbers
        for ii in range(len(columns)):
          try:
            lst = list(columns[ii].lower())
            while lst.count('d') > 0:
              lst[lst.index('d')] = 'e'
            columns[ii] = ''.join(lst)
            columns[ii] = float(columns[ii])
          except:
            raise Exception
        # Add to writing buffer
        outhe.append(line+'\n')
      # Write buffer
      for out in outhe:
          f.write(out)
    # Could not read the helium number density, so output none
    except:
      f.write('0\n')

  f.close()

if __name__ == "__main__":
    rAtmo()
