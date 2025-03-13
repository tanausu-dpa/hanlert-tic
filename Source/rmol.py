import sys, math, os

#####################
# rAtom()
#
# Tanaus\'u del Pino Alem\'an (IAC)
#
# 17/12/2024:  V4.0.0 - Changed global version (TdPA)
#
#####################

def rMol():
  ''' Reads the input molecule file specified as argument.
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
      print((msg+' in rmol.py'))
    else:
      # Check file exists
      exist = os.path.isfile(fil)
      # Open to write or append
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      # Write in file and close
      fv.write(msg+' in rmol.py\n')
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
    verbose(' # No molecule file found', verbfile, verbosity)
    filename = 'tmp_mol_'+dni
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
      lines_n.append(line)
  lines = lines_n

  # Start output file
  filename = 'tmp_mol_'+dni
  f = open(filename,'w')
  f.write('1\n')


  #
  # File build
  #

  # Label of the molecule
  iline = 0
  line = lines[iline]
  name = list(line.upper())
  f.write('{0}\n'.format(len(name)))
  line = ''.join(name)
  f.write(line)
  f.write('\n')

  #  Mass of the molecule
  iline += 1
  line = lines[iline]
  lst = list(line.lower())
  while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
  line = ''.join(lst)
  mass = float(line.strip())
  f.write(line)
  f.write('\n')

  #  Charge of the molecule
  iline += 1
  line = lines[iline]
  lst = list(line.lower())
  while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
  line = ''.join(lst)
  z = int(line.strip())
  f.write(line)
  f.write('\n')

  # List of atoms
  numbers = ['1','2','3','4','5','6','7','8','9']
  letters = ['A','B','C','D','E','F','G','H','I','J', \
             'K','L','M','N','O','P','Q','R','S','T', \
             'U','V','W','X','Y','Z']
  iline += 1
  line = lines[iline]
  atoms = line.strip().split(',')
  na = len(atoms)
  f.write('{0}\n'.format(len(atoms)))
  for atm in atoms:
    hadnum = False
    strlet = False
    num = ''
    name = ''
    pieces = list(atm)
    for pie in pieces:
      if pie in numbers:
        if strlet:
          verbose(' # Number after letter in elements of ' + \
                  'molecule, wrong format', verbfile, verbosity)
          abort(f,filename)
        hadnum = True
        num += pie
      if pie.upper() in letters:
        strlet = True
        if not hadnum:
          num += '1'
          hadnum = True
        name += pie
    if(len(list(name)) < 2):
      name = [' '] + [name[0]]
      name = name[0:2]
      name = ''.join(name)
    f.write('{0}\n'.format(num))
    f.write('{0}\n'.format(name.upper()))

  #  Dissociation energy
  iline += 1
  line = lines[iline]
  lst = list(line.lower())
  while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
  line = ''.join(lst)
  den = float(line.strip())
  f.write(line)
  f.write('\n')

  # Partition function type
  if z == 0:
    pft = ['KURUCZ_70','KURUCZ_85','SAUVAL_TATUM_84', \
           'IRWIN_81','TSUJI_73']
  elif z == 1:
    pft = ['KURUCZ_70','KURUCZ_85']
  else:
    verbose(' # Charge must be 0 or 1', verbfile, verbosity)
    abort(f,filename)
  dpft = {'KURUCZ_70':'0', 'KURUCZ_85':'1', \
          'SAUVAL_TATUM_84':'2','IRWIN_81':'3', \
          'TSUJI_73':'4'}
  iline += 1
  line = lines[iline].strip().upper()
  check = True
  for pf in pft:
    if line == pf:
      check = False
  if check:
    msg = ' # Cannot recognize partition function type, ' + \
          'for a molecule with charge {0} the '.format(z) + \
          'options are:'
    for pf in pft:
        msg += ' '+pf
    verbose(msg,verbfile,verbosity)
    abort(f,filename)
  f.write(dpft[line])
  f.write('\n')

  # Tmin and Tmax
  iline += 1
  line = lines[iline]
  tmps = line.split()
  if len(tmps) != 2:
    verbose(' # Wrong number of temperature inputs, ' + \
            'expected 2, got ' + \
            '{0}'.format(len(tmps)), verbfile, verbosity)
    abort(f,filename)
  for ii in range(len(tmps)):
    lst = list(tmps[ii].lower())
    while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
    tmps[ii] = ''.join(lst)
    tmps[ii] = float(tmps[ii])
  f.write(line)
  f.write('\n')

  # Partition function coefficients
  iline += 1
  line = lines[iline]
  pfc = line.split()
  npfc = int(pfc[0])
  pfc.pop(0)
  if npfc != len(pfc):
    verbose(' # Expected ' + \
            '{0}'.format(npfc) + \
            'coefficients, but only got ' + \
            '{0}'.format(len(pfc)),verbfile,verbosity)
    abort(f,filename)
  f.write('{0}\n'.format(npfc))
  line = ' '.join(pfc)
  if npfc > 0:
    f.write(line)
    f.write('\n')

  # Equilibrium coefficients
  iline += 1
  line = lines[iline]
  efc = line.split()
  nefc = int(efc[0])
  efc.pop(0)
  if nefc != len(efc):
    verbose(' # Expected ' + \
            '{0}'.format(nefc) + \
            'coefficients, but only got ' + \
            '{0}'.format(len(efc)),verbfile,verbosity)
    abort(f,filename)
  f.write('{0}\n'.format(nefc))
  line = ' '.join(efc)
  if nefc > 0:
    f.write(line)
    f.write('\n')

  # Close file
  f.close()

if __name__ == "__main__":
    rMol()

