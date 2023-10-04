import sys, math, os

#####################
# rAtmo()
#
# Tanaus\'u del Pino Alem\'an
# Ricky Egeland
#
# 08/07/2023:  V3.0.2 - Bugfix: Wrong indents (TdPA)
#
# 07/03/2023:  V3.0.1 - Now the type of model atmosphere is
#                       written in the temporal file before
#                       the thermal variables (TdPA)
#
# 06/29/2022:  V3.0.0 - Changed global version (TdPA)
#
# 03/17/2021:  V2.0.0 - Changed global version (TdPA)
#
# 01/13/2021 : V1.1.5 - Added proper abortion in abort() (TdPA)
#                     - Corrected typo in error message (TdPA)
#
# 06/05/2020 : V1.1.4 - Python 3 compatible (RE)
#
# 10/12/2019 : V1.1.3 - Bugfix: Improved verbosity of error for
#                       wrong reference wavelength (TdPA)
#
# 09/24/2019 : V1.1.2 - Added alternative inputs for the density
#                       specification (TdPA)
#
# 09/13/2019 : V1.1.1 - Added tau scale to the possible inputs (TdPA)
#                     - tau scales also admit specifying the reference
#                       wavelength (TdPA)
#
# 02/14/2019 : V1.1.0 - Improved verbosity (TdPA)
#
# 09/14/2017:  V1.0.1 - Added an ID to the files (TdPA)
#
# 04/17/2017:  V1.0.0 - First version (TdPA)
#
#####################

def rAtmo():
  ''' Reads the input atmospheric file specified as argument.
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
      print((msg+' in ratmo.py'))
    else:
      exist = os.path.isfile(fil)
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      fv.write(msg+' in ratmo.py\n')
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

  # Output file
  filename = 'tmp_atmo_'+dni
  f = open(filename,'w')
  f.write('1\n')

  # File build
  iline = 0
  line = lines[iline].strip()
  name = line

  # Type of scale
  iline += 1
  line = lines[iline].strip()
  scale = line.split()
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
  else:
    tf = '0.2\n'
  if 'HEIGHT' in scale[0]:
    f.write('H\n')
  elif 'TAU' in scale[0]:
    f.write('T\n')
 #elif 'MASS' in scale:
 #  f.write('M\n')
  else:
    verbose(' # Scale not recognized in model atmosphere', \
            verbfile, verbosity)
    abort(f,filename)
  f.write(tf)

  # Log(g), not used
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
  #f.write('{0:6.2f}\n'.format(logg))

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
  #f.write('{0:5d}'.format(NZ))

  # Initialize atmo to write
  thermo = []
  dens = []
  nd = False

  # Thermodynamical quantities
  for iz in range(NZ):
    iline += 1
    line = lines[iline].strip()
    columns = line.split()
    if len(columns) < 5:
      verbose(' # Needed at least five columns in ' + \
              'data block 1, got {0}'.format(len(colums)), \
              verbfile, verbosity)
      abort(f,filename)
    while len(columns) < 7:
      columns.append('0')
      line += ' 0'
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
    thermo.append(line+'\n')
    #form = '{22.16e}'*5 + '\n'
    #f.write(form.format(columns))

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
    f.write('0\n')
    nd = True

    # Go back one line
    iline -= 1

    # Hydrogen densities
    for iz in range(NZ):
      iline += 1
      line = lines[iline].strip()
      columns = line.split()
      if (len(columns)) < 2:
        verbose(' # Needed at least two columns in ' + \
                'data block 2, got {0}'.format(len(columns)), \
                verbfile, verbosity)
        abort(f,filename)
      changed = False
      while (len(columns)) < 6:
        columns.insert(-1,'0')
        changed = True
      if changed:
        line = ' '.join(columns)
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
      #form = '{22.16e}'*6 + '\n'
      #f.write(form.format(columns))

  # Output thermo
  for line in thermo:
    f.write(line)
  # Output number densities
  for line in dens:
    f.write(line)

  # If number densities
  if nd:

    # Helium densities
    try:
      outhe = ['1']
      for iz in range(NZ):
        iline += 1
        line = lines[iline].strip()
        columns = line.split()
        if (len(columns)) < 3:
          raise Exception
        changed = False
        while (len(columns)) < 4:
          columns.insert(1,'0')
          changed = True
        if changed:
          line = ' '.join(columns)
        for ii in range(len(columns)):
          try:
            lst = list(columns[ii].lower())
            while lst.count('d') > 0:
              lst[lst.index('d')] = 'e'
            columns[ii] = ''.join(lst)
            columns[ii] = float(columns[ii])
          except:
            raise Exception
        outhe.append(line+'\n')
      for out in outhe:
          f.write(out)
    except:
      f.write('0\n')

  f.close()

if __name__ == "__main__":
    rAtmo()
