import sys, math, os

#####################
# rAtom()
#
# Tanaus\'u del Pino Alem\'an
# Ricky Egeland
#
# 10/04/2023:  V3.0.1 - Just formatting, not advancing version (TdPA)
#
# 04/25/2023:  V3.0.1 - Bugfix: The float check in the collisions
#                       was happening before fully correcting for
#                       python suitable format (TdPA)
#
# 06/29/2022:  V3.0.0 - Changed global version (TdPA)
#                     - Removed the possibility of explicit
#                       input for depolarizing collisions (TdPA)
#
# 09/30/2021:  V2.0.1 - Prepare to deal with forbidden lines (TdPA)
#
# 03/17/2021:  V2.0.0 - Changed global version (TdPA)
#
# 01/13/2021 : V1.3.6 - Added proper abortion in abort() (TdPA)
#
# 09/28/2020 : V1.3.5 - Change defaults of unsold (TdPA)
#
# 06/05/2020 : V1.3.4 - Python 3 compatible (RE)
#
# 08/14/2019 : V1.3.3 - Bugfix: Indexes in coldic dictionary had
#                       to be strings, not integers (TdPA)
#
# 07/23/2019 : V1.3.2 - More possible flags for the type of
#                       photoionization: neutral, ion, lneutral,
#                       lion, sneutral, sion, none, linear, and
#                       spline (TdPA)
#
# 03/06/2019 : V1.3.1 - Now there is an extra integer in the radiative
#                       transition line in atomic models, and it is
#                       optional (TdPA)
#
# 02/14/2019 : V1.3.0 - Improved verbosity (TdPA)
#
# 11/02/2017:  V1.2.1 - Multilevel expects a Land\'e factor
#
# 10/30/2017:  V1.2.0 - Changed how the multilevel works
#
# 10/11/2017:  V1.1.0 - Changed the structure of the energy part
#
# 09/22/2017:  V1.0.3 - Removed some nested index duplicity
#
# 09/14/2017:  V1.0.2 - Added an ID to the files
#
# 07/19/2017:  V1.0.1 - Typo in one error mesagge
#                     - No message for No collisions (moved to main
#                       code)
#
# 04/18/2017:  V1.0.0 - First version
#
#####################

def rAtom():
  ''' Reads the input atomic file specified as argument.
  '''

  def abort(f,name):
    f = open(name,'w')
    f.write('-1')
    f.close()
    sys.exit()

  def verbose(msg, fil, verb):

    # If being verbose
    if (verb):
      print((msg+' in ratom.py'))
    else:
      exist = os.path.isfile(fil)
      if (exist):
        fv = open(fil,'a')
      else:
        fv = open(fil,'w')
      fv.write(msg+' in ratom.py\n')
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
    verbose(' # No atomic file found', verbfile, verbosity)
    filename = 'tmp_atom_'+dni
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

  # Output file
  filename = 'tmp_atom_'+dni
  output = ['1\n']

  # File build

  #  Label of the atom
  iline = 0
  line = lines[iline]
  name = list(line.upper())
  if(len(name) < 2):
    name = [' '] + [name[0]]
  name = name[0:2]
  line = ''.join(name)
  output.append(line)
  output.append('\n')

  #  Mass of the atom in AMU
  iline += 1
  line = lines[iline]
  lst = list(line.lower())
  while lst.count('d') > 0:
    lst[lst.index('d')] = 'e'
  line = ''.join(lst)
  rmass = float(line.strip())
  output.append(line)
  output.append('\n')

  #  Abundance, isotopic % and renormalization
  iline += 1
  line = lines[iline]
  cols = line.split()
  if len(cols) < 3:
    cols.append('1')
  for ii in range(len(cols)):
    col = cols[ii]
    lst = list(col.lower())
    while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
    col = ''.join(lst)
    if ii < 2:
      col = float(col)
    else:
      col = int(col)
    output.append('{0}'.format(col))
    output.append('\n')

  #  Type of model
  models = ['MULTILEVEL','MULTITERM','MLEVEL','MTERM','ML','MT']
  iline += 1
  line = lines[iline]
  if line.upper().strip() not in models:
    verbose(' # Model type not recognized. Acceted models ' + \
            'are {0} and {1}'.format(*models[0:1]), verbfile, \
            verbosity)
    abort(f,filename)
  if line.upper().strip() in models[0] or \
     line.upper().strip() in models[2] or \
     line.upper().strip() in models[4]:
    model = 1
  else:
    model = 0
  output.append('{0}'.format(model))
  output.append('\n')

  #  Number of levels, lines, photoionizations and depol. col.
  iline += 1
  line = lines[iline]
  dims = line.split()
  if len(dims) != 4:
    verbose(' # Wrong number of dimension, ' + \
            'expected 4, got {0}'.format(len(dims)), \
            verbfile,verbosity)
    abort(f,filename)
  for ii in range(len(dims)):
    lst = list(dims[ii].lower())
    while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
    dims[ii] = ''.join(lst)
    dims[ii] = int(dims[ii])
  if len(dims) <= 2:
    verbose(' # You need at least 2 terms', verbfile, verbosity)
    abort(f,filename)
  if len(dims) <= 1:
    verbose(' # You need at least 1 transition', verbfile, \
            verbosity)
    abort(f,filename)
  idim = len(output)
  output.append(line)
  output.append('\n')

  # Find maximum nJ and write it
  if model == 1:
    output.append('1\n')
  else:
    iline_0 = iline
    nJmax = 0
    for ii in range(dims[0]):
      iline += 1
      line = lines[iline]
      q_dims = line.split()
      if len(q_dims) != 4:
        verbose(' # Wrong number of term inputs, ' + \
                'expected 4, got {0} for term {1}'.format( \
                len(q_dims),ii+1), verbfile, verbosity)
        abort(f,filename)
      for kk in range(len(q_dims)):
        lst = list(q_dims[kk].lower())
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        q_dims[kk] = ''.join(lst)
        q_dims[kk] = float(q_dims[kk])
      rJmin = abs(q_dims[0] - q_dims[1])
      rJmax = q_dims[0] + q_dims[1]
      nJ = int(rJmax - rJmin) + 1
      if nJ > nJmax:
        nJmax = nJ
      for jj in range(nJ):
        iline += 1
    output.append(str(nJmax)+'\n')
    iline = iline_0

  # Energy levels
  nJv = []
  # Multiterm
  if model == 0:
    for ii in range(dims[0]):
      iline += 1
      line = lines[iline]
      q_dims = line.split()
      if len(q_dims) != 4:
        verbose(' # Wrong number of term inputs, ' + \
                'expected 4, got {0} for term {1}'.format( \
                len(q_dims),ii+1), verbfile, verbosity)
        abort(f,filename)
      for kk in range(len(q_dims)):
        lst = list(q_dims[kk].lower())
        while lst.count('d') > 0:
            lst[lst.index('d')] = 'e'
        q_dims[kk] = ''.join(lst)
        q_dims[kk] = float(q_dims[kk])
      output.append(line)
      output.append('\n')
      rJmin = abs(q_dims[0] - q_dims[1])
      rJmax = q_dims[0] + q_dims[1]
      nJ = int(rJmax - rJmin) + 1
      nJv.append(nJ)
      for jj in range(nJ):
        iline += 1
        line = lines[iline]
        eners = line.split()
        if len(eners) != 4:
          verbose(' # Wrong number of level inputs, ' + \
                  'expected 4, got ' + \
                  '{0} for level {1} of term {2}'.format( \
                  len(eners),jj+1,ii+1), verbfile, verbosity)
          abort(f,filename)
        for kk in range(len(eners)):
          lst = list(eners[kk].lower())
          while lst.count('d') > 0:
            lst[lst.index('d')] = 'e'
          eners[kk] = ''.join(lst)
          if kk < 3:
            eners[kk] = float(eners[kk])
          else:
            eners[kk] = int(eners[kk])
        output.append(line)
        output.append('\n')
  # Multilevel
  else:
    for ii in range(dims[0]):
      iline += 1
      line = lines[iline]
      q_dims = line.split()
      if len(q_dims) != 4:
        verbose(' # Wrong number of term inputs, ' + \
                'expected 4, got {0} for term {1}'.format( \
                len(q_dims),ii+1), verbfile, verbosity)
        abort(f,filename)
      for kk in range(len(q_dims)):
        lst = list(q_dims[kk].lower())
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        q_dims[kk] = ''.join(lst)
        q_dims[kk] = float(q_dims[kk])
      output.append(line)
      output.append('\n')
      nJv.append(1)
      iline += 1
      line = lines[iline]
      eners = line.split()
      if len(eners) != 5:
        verbose(' # Wrong number of level inputs, ' + \
                'expected 5, got ' + \
                '{0} for level {1}'.format( \
                len(eners),ii+1), verbfile, verbosity)
        abort(f,filename)
      for kk in range(len(eners)):
        lst = list(eners[kk].lower())
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        eners[kk] = ''.join(lst)
        if kk < 4:
          eners[kk] = float(eners[kk])
        else:
          eners[kk] = int(eners[kk])
      output.append(line)
      output.append('\n')

  # Transition lines
  bark = ['s','p','d','f']
  bardic = {'s':'0','p':'1','d':'2','f':'3'}
  broads = ['barklem','unsold','param']
  broads_dic = {'barklem':'0','unsold':'1','param':'2'}
  classe = ['e1','m1','e2','m2','un']
  classe_dic = {'e1':1,'m1':2,'e2':3,'m2':4,'un':5}
  nfor = 0
  for ii in range(dims[1]):
    iline += 1
    line = lines[iline]
    cols = line.lower().split()
    # Check format
    if cols[3] in classe:
      new = True
      numbers = [0,1,2,6,8,9,10,11,12,13,14,15]
    elif cols[3] in broads:
      new = False
      numbers = [0,1,2,5,7,8,9,10,11,12,13,14]
    else:
      verbose(' # The fourth element in transition ' + \
            '{0}'.format(ii+1) + \
            ', {0}, is not a broadening of the '.format(cols[3]) + \
            'following: {0}, {1}, or {2},'.format(*broads) + \
            'nor a transition type of the '.format(cols[3]) + \
            'following: {0}, {1}, {2}, {3} or {4},'.format(*classe), \
              verbfile, verbosity)
      abort(f,filename)
    if new:
      if len(cols) < 15:
        verbose(' # Wrong number of transition inputs, ' + \
                'expected at least 14, got ' + \
                '{0} for transition {1}'.format(len(cols),ii+1), \
                verbfile, verbosity)
        abort(f,filename)
      if len(cols) > 16:
        verbose(' # Wrong number of transition inputs, ' + \
                'expected not more than 15, got ' + \
                '{0} for transition {1}'.format(len(cols),ii+1), \
                verbfile, verbosity)
        abort(f,filename)
      if len(cols) != 16:
          cols.append('1')
      if cols[3].lower() not in classe:
        verbose(' # The transition type in transition ' + \
            '{0}'.format(ii+1) + \
            ', {0}, is none of the '.format(cols[3]) + \
            'following: {0}, {1}, {2}, {3}, or {4}'.format(*classe), \
                verbfile, verbosity)
        abort(f,filename)
      if cols[3].lower() != classe[0] and \
         cols[3].lower() != classe[1]:
          nfor += 1
      if cols[4] not in broads:
        verbose(' # The broadening type in transition ' + \
                '{0}'.format(ii+1) + \
                ', {0}, is none of the '.format(cols[4]) + \
                'following: {0}, {1}, or {2}'.format(*broads), \
                verbfile, verbosity)
        abort(f,filename)
      if cols[4] == 'barklem':
        changed = False
        if cols[5] not in bark:
          verbose(' # There is no Barklem data for l: ' + \
                  '{0}'.format(cols[5]) + \
                  ', changed to Unsold without enhancement for' + \
                  ' line {0}'.format(ii+1), verbfile, verbosity)
          cols[4] = 'unsold'
          cols[5] = '1e0'
          cols[6] = '0e0'
          cols[7] = '1e0'
          cols[8] = '0e0'
          changed = True
        if cols[7] not in bark and not changed:
          verbose(' # There is no Barklem data for l: ' + \
                  '{0}'.format(cols[7]) + \
                  ', changed to Unsold without enhancement for' + \
                  ' line {0}'.format(ii+1), verbfile, verbosity)
          cols[4] = 'unsold'
          cols[5] = '1e0'
          cols[6] = '0e0'
          cols[7] = '1e0'
          cols[8] = '0e0'
        if not changed:
          cols[5] = bardic[cols[5]]
          cols[7] = bardic[cols[7]]
      if cols[4] == 'unsold' or cols[4] == 'param':
        numbers.append(5)
        numbers.append(7)
      for kk in numbers:
        lst = list(cols[kk])
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        cols[kk] = ''.join(lst)
        num = float(cols[kk])
      output.append(' '.join(cols[0:3]))
      output.append('\n')
      output.append(classe_dic[cols[3]])
      output.append('\n')
      output.append(broads_dic[cols[4]])
      output.append('\n')
      for jj in range(5,9):
        output.append(cols[jj] + ' ')
        output.append('\n')
      output.append(' '.join(cols[9:]))
      output.append('\n')
    # Old
    else:
      if len(cols) < 14:
        verbose(' # Wrong number of transition inputs, ' + \
                'expected at least 14, got ' + \
                '{0} for transition {1}'.format(len(cols),ii+1), \
                verbfile, verbosity)
        abort(f,filename)
      if len(cols) > 15:
        verbose(' # Wrong number of transition inputs, ' + \
                'expected not more than 15, got ' + \
                '{0} for transition {1}'.format(len(cols),ii+1), \
                verbfile, verbosity)
        abort(f,filename)
      if len(cols) != 15:
          cols.append('1')
      if cols[3] not in broads:
        verbose(' # The broadening type in transition ' + \
                '{0}'.format(ii+1) + \
                ', {0}, is none of the '.format(cols[3]) + \
                'following: {0}, {1}, or {2}'.format(*broads), \
                verbfile, verbosity)
        abort(f,filename)
      if cols[3] == 'barklem':
        changed = False
        if cols[4] not in bark:
          verbose(' # There is no Barklem data for l: ' + \
                  '{0}'.format(cols[4]) + \
                  ', changed to Unsold without enhancement for' + \
                  ' line {0}'.format(ii+1), verbfile, verbosity)
          cols[3] = 'unsold'
          cols[4] = '1e0'
          cols[5] = '0e0'
          cols[6] = '1e0'
          cols[7] = '0e0'
          changed = True
        if cols[6] not in bark and not changed:
          verbose(' # There is no Barklem data for l: ' + \
                  '{0}'.format(cols[6]) + \
                  ', changed to Unsold without enhancement for' + \
                  ' line {0}'.format(ii+1), verbfile, verbosity)
          cols[3] = 'unsold'
          cols[4] = '1e0'
          cols[5] = '0e0'
          cols[6] = '1e0'
          cols[7] = '0e0'
        if not changed:
          cols[4] = bardic[cols[4]]
          cols[6] = bardic[cols[6]]
      if cols[3] == 'unsold' or cols[3] == 'param':
        numbers.append(4)
        numbers.append(6)
      for kk in numbers:
        lst = list(cols[kk])
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        cols[kk] = ''.join(lst)
        num = float(cols[kk])
      output.append(' '.join(cols[0:3]))
      output.append('\n')
#     output.append('1\n')
      output.append(broads_dic[cols[3]])
      output.append('\n')
      for jj in range(4,8):
        output.append(cols[jj] + ' ')
        output.append('\n')
      output.append(' '.join(cols[8:]))
      output.append('\n')
# output[idim] += ' {0:14d}'.format(nfor)

  # Depolarizing collisions
  depcol = ['fit']
  for ii in range(dims[3]):
    iline += 1
    line = lines[iline]
    cols = line.lower().split()
    if len(cols) != 2:
      verbose(' # Wrong number of depol. inputs, ' + \
              'expected 2 indicators, got ' + \
              '{0} for entry {1}'.format(len(cols),ii+1), \
              verbfile, verbosity)
      abort(f,filename)
    for col in cols:
      nentry = int(col)
    output.append(line)
    output.append('\n')
    for jj in range(nentry):
      iline += 1
      line = lines[iline]
      cols = line.lower().split()
      if cols[1] not in depcol:
        verbose(' # The input type in depol. col. ' + \
                '{0} of input {1}'.format(jj+1,ii+1) + \
                ', {0}, is none of the '.format(cols[1]) + \
                'following: {0}, or {1}'.format(*depcol), \
                verbfile, verbosity)
        abort(f,filename)
      num = int(cols[0])
      output.append(cols[0])
      output.append('\n')
      if cols[1] == 'fit':
        output.append('f')
      output.append('\n')
      if cols[1] == 'fit':
        output.append('1')
        output.append('\n')
        iline += 1
        line = lines[iline]
        cols = line.lower().split()
        if len(cols) < 2:
          verbose(' # Wrong number of depol. inputs, ' + \
                  'expected at least 2 parameters, got ' + \
                  '{0} for entry '.format(len(cols)) + \
                  '{0} of input {1}'.format(jj+1,ii+1), \
                  verbfile, verbosity)
          abort(f,filename)
        for col in cols:
          lst = list(col)
          while lst.count('d') > 0:
            lst[lst.index('d')] = 'e'
          col = ''.join(lst)
        output.append(' '.join(cols))
        output.append('\n')

  # Photoionization cross section
  photo = ['hydrogenic','explicit']
  for ii in range(dims[2]):
    iline += 1
    line = lines[iline]
    cols = line.lower().split()
    if len(cols) != 4:
      verbose(' # Wrong number photoionization inputs, ' + \
              'expected 4 indicators, got ' + \
              '{0} for entry {1}'.format(len(cols),ii+1), \
              verbfile, verbosity)
      abort(f,filename)
    if cols[2] not in photo:
      verbose(' # The input type in photoionization ' + \
              'entry {0}'.format(ii+1) + \
              ', {0}, is none of the '.format(cols[2]) + \
              'following: {0}, or {1}'.format(*photo), \
              verbfile, verbosity)
      abort(f,filename)
    num = int(cols[0])
    num = int(cols[1])
    output.append(' '.join(cols[0:2]))
    output.append('\n')
    num = int(cols[3])
    if cols[2] == 'explicit':
      nentry = num
      output.append('e')
    if cols[2] == 'hydrogenic':
      nentry = 1
      output.append('h')
    output.append('\n')
    output.append(cols[3])
    output.append('\n')
    for jj in range(nentry):
      iline += 1
      line = lines[iline]
      cols = line.lower().split()
      if len(cols) != 2:
        verbose(' # Wrong number photoionization data, ' + \
                'expected 2 numbers, got ' + \
                '{0} for entry {1}'.format(len(cols),jj+1) + \
                ' of input {0}'.format(ii+1), verbfile, verbosity)
        abort(f,filename)
      for col in cols:
        lst = list(col)
        while lst.count('d') > 0:
          lst[lst.index('d')] = 'e'
        col = ''.join(lst)
      output.append(' '.join(cols))
      output.append('\n')

  # Collisional rates
  try:
    rt = False
    coldic = {'be':'0', 'fe':'1', 'bp':'2', 'bh':'3', 'fh':'4', \
              'c0':'5', 'c+':'6'}
    ctype = ['be','fe','bp','bh','fh','c0','c+']
    nlins = len(lines) - iline - 1
    ncol = 0
    nclin = 0
    outcol = []
    for ii in range(nlins):
      iline +=1
      line = lines[iline]
      cols = line.lower().split()
      flag1 = cols[0]
      if(flag1.lower() == 'cols'):

        NT = len(cols) - 3
        if NT < 2:
          verbose(' # Wrong number collisions data, ' + \
                  'expected at least 5 inputs, got ' + \
                  '{0}'.format(len(cols)), verbfile, verbosity)
          abort(f,filename)
        flag1 = -1
        flag2 = cols[1]
        flag3 = cols[2]

        # Check flag2
        if(flag2.lower() in 'term'):
          flag2 = 0
        elif(flag2.lower() in 'level'):
          flag2 = 1
        else:
          verbose(' # Wrong flag2 in collisions', verbfile, \
                  verbosity)
          flag2 = -1
          raise Exception

        # Check flag 3
        if(flag3.lower() in 'neutral'):
          flag3 = 0
        elif(flag3.lower() in 'ion'):
          flag3 = 1
        elif(flag3.lower() in 'sneutral'):
          flag3 = 0
        elif(flag3.lower() in 'sion'):
          flag3 = 1
        elif(flag3.lower() in 'lneutral'):
          flag3 = 2
        elif(flag3.lower() in 'lion'):
          flag3 = 3
        elif(flag3.lower() in 'linear'):
          flag3 = -2
        elif(flag3.lower() in 'spline'):
          flag3 = -1
        elif(flag3.lower() in 'none'):
          flag3 = -1
        else:
          msg = ' # Warning: unknown flag3 {0} in collisions'
          msg = msg.format(flag3)
          verbose(msg, verbfile, verbosity)
          flag3 = -1
          raise Exception

        # Get temperatures
        temps = cols[3:]
        outcol.append('{0}\n{1}\n{2}\n{3}'. \
                      format(flag1,flag2,flag3,NT))
        outcol.append(' '.join(temps))
        rt = True
        nclin += 1

      elif flag1.lower() in ctype :

        if not rt:
          verbose(' # No temperature axis defined before ' + \
                  'the collisions.', verbfile, verbosity)
          raise Exception

        outcol.append(coldic[flag1])

        if cols[-1] == 'f':
          outcol.append('1')
          cols.pop(-1)
        else:
          outcol.append('0')

        if flag2 == 0 and coldic[flag1] == '1':
          verbose(' # Ionizing collisions only admits ' + \
                  'level wise indexing', verbfile, verbosity)
          raise Exception


        if len(cols) != (3+NT):
          verbose(' # Wrong number collisional data, ' + \
                  'expected {0} '.format(NT+3) + \
                  ', got ' + \
                  '{0} for entry {1}.'.format(len(cols),ii+1), \
                  verbfile, verbosity)
          raise Exception

        num = int(cols[1])
        num = int(cols[2])
        outcol.append(' '.join(cols[1:3]))
        data = cols[3:]
        for dat in data:
          lst = list(dat)
          while lst.count('d') > 0:
              lst[lst.index('d')] = 'e'
          dat = ''.join(lst)
          flt = float(dat)
        outcol.append(' '.join(data))
        ncol += 1

      else:
        if cols[0] not in ctype:
            verbose(' # The flag in collisions ' + \
                    'entry {0}'.format(ii+1) + \
                    ', {0}, is none of the '.format(cols[0]) + \
                    'following: cols, {0}, or {1}'.format(*ctype), \
                    verbfile, verbosity)
            raise Exception

    output.append('{0}\n'.format(ncol))
    output.append('{0}\n'.format(ncol+nclin))
    for olin in outcol:
      output.append(olin)
      output.append('\n')
    if not rt:
      output.append('0\n')
  except ValueError:
    msg= ' # There was a value error when trying to write' + \
         ' the collisions. Maybe you have spurious characters'
    verbose(msg.format(), verbfile, verbosity)
    output.append('0\n')
  except:
    output.append('0\n')

  # Write file
  f = open(filename,'w')
  for ou in output:
    f.write(ou)
  f.close()

if __name__ == "__main__":
    rAtom()
