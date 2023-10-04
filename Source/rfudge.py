import sys, math, os

#####################
# rfudge()
#
# Tanaus\'u del Pino Alem\'an
# Ricky Egeland
#
# 06/29/2022:  V3.0.0 - Changed global version (TdPA)
#
# 03/17/2021:  V2.0.0 - Changed global version (TdPA)
#
# 01/13/2021 : V1.1.5 - Added proper abortion in abort() (TdPA)
#
# 06/05/2020 : V1.1.4 - Python 3 compatible (RE)
#
# 02/15/2019 : V1.1.0 - Improved verbosity (TdPA)
#
# 09/14/2017:  V1.0.1 - Added an ID to the files
#
# 04/19/2017:  V1.0.0 - First Version
#
#####################

def rfudge():
  ''' Reads the fudge file specified as argument.
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
      print((msg+' in rfudge.py'))
    else:
      exist = os.path.isfile(fil)
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      fv.write(msg+' in rfudge.py\n')
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
    verbose(' # No fudge file found', verbfile, verbosity)
    filename = 'tmp_fud_'+dni
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
  filename = 'tmp_fud_'+dni
  f = open(filename,'w')
  f.write('1\n')

  # Output buff
  output = []

  # Direction checker
  freq = []

  # File build
  for line in lines:
      cols = line.split()
      ncol = len(cols)
      if ncol < 2:
          continue
      brek = False
      for ele in cols:
          try:
              dump = float(ele)
          except:
              brek = True
              continue
      if brek:
          continue
      while ncol < 4:
          cols.append('1.00')
          ncol = len(cols)
      cols[0] = '{0}'.format(float(cols[0]))
      freq.append(cols[0])
      output.append(' '.join(cols[:5]))

  if float(freq[0]) > float(freq[1]):
    output = output[::-1]

  f.write('{0}\n'.format(len(output)))
  for out in output:
    f.write(out+'\n')

  f.close()

if __name__ == "__main__":
    rfudge()
