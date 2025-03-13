import sys, math, os

#####################
# rfudge()
#
# Tanaus\'u del Pino Alem\'an (IAC)
#
# 13/12/2024:  V4.0.0 - Changed global version (TdPA)
#
#####################

def rfudge():
  ''' Reads the fudge file specified as argument.
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
      print((msg+' in rfudge.py'))
    else:
      # Check file exists
      exist = os.path.isfile(fil)
      # Open to write or append
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      # Write in file and close
      fv.write(msg+' in rfudge.py\n')
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

  # Start output file
  filename = 'tmp_fud_'+dni
  f = open(filename,'w')
  f.write('1\n')

  # Output buffer
  output = []

  # Direction checker
  freq = []

  #
  # File build
  #

  # For each read line
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

  # If reverse wavelength order
  if float(freq[0]) > float(freq[1]):
    output = output[::-1]

  # Output in file
  f.write('{0}\n'.format(len(output)))
  for out in output:
    f.write(out+'\n')

  # Close file
  f.close()

if __name__ == "__main__":
    rfudge()
