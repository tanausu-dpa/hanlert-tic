import sys, math, os

#####################
# rBField()
#
# Tanaus\'u del Pino Alem\'an
# Ricky Egeland
#
# 06/29/2022:  V3.0.0 - Changed global version (TdPA)
#
# 03/17/2021:  V2.0.0 - Changed global version (TdPA)
#
# 01/13/2021 : V1.1.2 - Added proper abortion in abort() (TdPA)
#
# 06/05/2020 : V1.1.1 - Python 3 compatible (RE)
#
# 02/14/2019 : V1.1.0 - Improved verbosity (TdPA)
#
# 09/14/2017(US):  V1.0.1 - Added an ID to the files
#
# 04/17/2017(US):  V1.0.0 - First version
#
#####################


def rBField():
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
      print((msg+' in rbfield.py'))
    else:
      exist = os.path.isfile(fil)
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      fv.write(msg+' in rbfield.py\n')
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

  # Output file
  filename='tmp_bfield_'+dni
  f=open(filename,'w')
  f.write('1\n')

  # File build
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
  for i in range(NZ):
    val = lines[2+i].split()
    if(len(val) != 3):
      verbose(' # The format must have three columns', \
              verbfile, verbosity)
      abort(f, filename)
    for j in range(3):
      f.write('{0:22.16e} '.format(float(val[j])))
    f.write('\n')
  f.close()

if __name__ == "__main__":
    rBField()
