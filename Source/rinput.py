import sys, math, os, shutil

#####################
# rInput()
#
# Tanaus\'u del Pino Alem\'an (IAC)
# Hao Li (IAC/NSSCC)
#
# 30/04/2025:  V4.0.4 - Bugfix: BFIELD_INPUT was ignored in the
#                       inversion mode (TdPA)
#
#####################

# Transform fortran double "d" to python float "e"
def interpret(val):
    ''' Interpret double precision nomenclature
    '''
    lst = list(val.lower())
    while lst.count('d') > 0:
      lst[lst.index('d')] = 'e'
    return ''.join(lst)

# Convert list of float string into ordered list of floats
def order_float(vec):
    ''' Convert a list of strings into an ordered list of floats.
        Values must be unique
    '''
    try:
      array = []
      for val in vec:
        array.append(float(val))
      array.sort()
      for i in range(1,len(array)):
        if array[i] <= array[i-1]:
          return []
      return array
    except:
      return []

# Order 2D ranges by its first dimension avoiding duplicates
def process_pixels(pixels):
    ''' Order pixels by first dimension and avoid duplicates
    '''
    output = []
    aux = {}
    lst = []
    for pixel in pixels:
      if pixel[0] in aux:
        if pixel[1] not in aux[pixel[0]]:
          aux[pixel[0]].append(pixel[1])
      else:
        aux[pixel[0]] = [pixel[1]]
        lst.append(pixel[0])
    lst.sort()
    for ix in lst:
        for iy in aux[ix]:
            output.append([ix,iy])
    return output

# Process a list of ranges ensuring no overlapping
def unique_ranges(NL,ilow,iup,iff,tran):
    ''' Check that a list of ranges does not have intersections
    '''
    # Create an auxiliary dictionary and list
    aux = {}
    keys = []

    # Fill starting wavelength and order
    for i in range(NL):
        if iup[i] < tran[0] or ilow[i] > tran[1]: continue
        if ilow[i] < tran[0]: ilow[i] = tran[0]
        if iup[i] > tran[1]: iup[i] = tran[1]
        aux[ilow[i]] = i
        keys.append(ilow[i])

    # Not valid limits
    if len(keys) < 1:
        return True,0,[],[],[]

    # Order the keys
    keys.sort()

    # Create new outputs
    low = []
    up = []
    ff = []

    # For each key, in order, add the info
    for key in keys:
        ind = aux[key]
        low.append(ilow[ind])
        up.append(iup[ind])
        ff.append(iff[ind])

    # Check no intersections
    # For each interval
    for iran in range(NL-1):

        # Check if next range shares limits
        if (low[iran+1] >= low[iran] and low[iran+1] <= up[iran]) or \
           (up[iran+1] >= low[iran] and up[iran+1] <= up[iran]):

            # Problem
            return False,0,[],[],[]

    # Return ordered ranges
    return True, len(keys), low, up, ff

# Order a list of wavelength ranges and combine if they overlap
def Worder(NL,doublets):
    ''' Order and combine wavelength ranges
    '''

    # Create an auxiliary dictionary and list
    aux = {}
    keys = []

    # Create output
    ndoublets = []

    # Fill starting wavelength and order
    for idou,dou in zip(range(NL),doublets):
        aux[dou[0]] = idou
        keys.append(dou[0])

    # Order the keys
    keys.sort()

    # For each key, in order, add the range
    for key in keys:
        ndoublets.append(doublets[aux[key]])

    # Initialize changes
    change = True

    # Do while there are changes
    while change:

        # Just started, so no changes
        change = False

        # For each interval
        for iran in range(len(ndoublets)-1):

            # Check if next range shares limits
            if (ndoublets[iran+1][0] >= ndoublets[iran][0] and \
                ndoublets[iran+1][0] <= ndoublets[iran][1]) or \
               (ndoublets[iran+1][1] >= ndoublets[iran][0] and \
                ndoublets[iran+1][1] <= ndoublets[iran][1]):

                # We are going to change the doublets
                change = True

                # We will also be removing one element
                NL -= 1

                # Combine in iran + 1 both ranges
                ndoublets[iran+1][0] = ndoublets[iran][0]
                ndoublets[iran+1][1] = max([ndoublets[iran][1], \
                                            ndoublets[iran+1][1]])
                # Drop iran
                ndoublets.pop(iran)

                # And start again
                break

    # Return true ranges
    return NL, ndoublets

def rInput():
  ''' Reads the input file specified as argument. It follows the
      format 'Keyword = Value', with the valid keywords hardcoded
      in this routine
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
  def verbose(msg, folder, verb):
    # If being verbose
    if (verb):
      # Just print
      print((msg+' in rinput.py'))
    else:
      # Check file exists
      exist = os.path.isfile(folder+'/verbose')
      # Open to write or append
      if (exist):
        fv = open(folder+'/verbose','a')
      else:
        fv = open(folder+'/verbose','w')
      # Write in file and close
      fv.write(msg+' in rinput.py\n')
      fv.close()

  # Routine to process an LTE line input in Kurucz format
  def process_LTEline_entry_kur(entry):
    ''' Process an entry for LTE lines in Kurucz format
    '''

    # Check exact size
    if len(entry) != 160:

        # Check first
        cols = entry.split()
        if len(cols[0]) != 11:
            entry = ' '*(11 - len(cols[0])) + entry

        # Check again
        if len(entry) != 160:

            # Add iso
            entry = entry + '------'

    # Last check
    if len(entry) != 160:
        return False,[]

    # Try processing
    try:

      # Initialize output
      lout = [1]

      # Atom and stage
      code = entry[18:24]
      atom, stage = code.split('.')
      atom = int(atom)
      stage = int(stage) + 1
      lout.append(atom)
      lout.append(stage)

      # Atomic quantities
      E1 = float(entry[24:36])*1e-5
      J1 = float(entry[36:41])
      g1 = float(entry[144:149])*1e-3
      E2 = float(entry[52:64])*1e-5
      J2 = float(entry[64:69])
      g2 = float(entry[149:154])*1e-3

      # E1 lower
      if E1 < E2:
          lout.append(E1)
          lout.append(J1)
          lout.append(g1)
          lout.append(E2)
          lout.append(J2)
          lout.append(g2)
      # E1 larger
      else:
          lout.append(E2)
          lout.append(J2)
          lout.append(g2)
          lout.append(E1)
          lout.append(J1)
          lout.append(g1)

      # Einstein
      code = entry[11:18] # loggf
      Aul = 10e0**float(code)
      Aul = Aul/(2e0*lout[-2]+1e0)
      # Constants
      e0 = 1.60217646e-19      #C(A*s) Fundamental charge in SI
      ep0 = 8.854187817e-12    #F*m**-1
      me = 9.10938188e-31      #Kg electron mass
      cl = 299792458e0         #m/s speed of light
      lamb = 1e2/(lout[-3] - lout[-6])
      Aul = 2e0*math.pi*e0*e0*Aul*1e10/(lamb*lamb*cl*ep0*me)
      lout.append(Aul)

      # VdW
      lout.append('3')
      code = entry[92:98]
      num = 10e0**float(code)
      lout.append(num)
      lout.append(0.)
      lout.append(0.)
      lout.append(0.)

      # Stark
      code = entry[86:92]
      num = -10e0**float(code)
      lout.append(num)

      # Collisional
      lout.append(0e0)

      # Radiative
      code = entry[80:86]
      num = 10e0**float(code)
      lout.append(num)

      # Frequencies
      lout.append(41)
      lout.append(31)
      lout.append(10.)
      lout.append(5.)
      lout.append('1')
      lout.append('-1e0')
      lout.append('1e90')

      # If we are here, everything is valid
      return True, lout

    except:
      return False, []

  # Routine to process an LTE line input in HanleRT format
  def process_LTEline_entry(entry,ofolder,verbosity):
    ''' Process an entry for LTE lines
    '''

    # Error
    def lerror(msg,ofolder,verbosity):
      verbose(' # LTE_LINE entry invalid ' + msg + ' ignored', \
              ofolder, verbosity)

    # Split in columns
    cols = entry.upper().split()

    # Check number
    if len(cols) < 21:

      # Not valid entry
      lerror('number of colums (<21)',ofolder,verbosity)
      return False, []

    # No Tlim or taulim
    elif len(cols) < 22:

      # Add fake limits
      cols += ['e1','-1e0','1e90']

    # No Tlim or taulim
    elif len(cols) < 23:

      # Add fake limits
      cols += ['-1e0','1e90']

    # No Tlim
    elif len(cols) < 24:

      # Add fake limit
      cols.append('1e90')

    #
    # Atom
    #
    i = 0
    msg = 'atomic label '+cols[i]
    try:
      cols[i] = int(cols[i])
      if cols[i] < 1 or cols[i] > 99:
        lerror(msg,ofolder,verbosity)
      iatom_type = 1
    except ValueError:
      if len(cols[i]) > 3:
        lerror(msg,ofolder,verbosity)
        return False, []
      iatom_type = 0
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    #
    # Stage
    #
    i = 1
    msg = 'stage '+cols[i]
    sdic = {'I': 1, \
            'II': 2, \
            'III': 3, \
            'IV': 4, \
            'V': 5, \
            'VI': 6, \
            'VII': 7, \
            'VIII': 8, \
            'IX': 9, \
            'X': 10, \
            'XI': 11, \
            'XII': 12, \
            'XIII': 13, \
            'XIV': 14, \
            'XV': 15, \
            'XVI': 16, \
            'XVII': 17, \
            'XVIII': 18, \
            'XIX': 19, \
            'XX': 20, \
            'XXI': 21, \
            'XXII': 22, \
            'XXIII': 23, \
            'XXIV': 24, \
            'XXV': 25}
    try:
      cols[i] = int(cols[i])
      if cols[i] < 0 or cols[i] > 99:
        lerror(msg,ofolder,verbosity)
        return False, []
      if iatom_type == 1:
        if cols[i] > cols[i-1]:
          lerror(msg,ofolder,verbosity)
          return False, []
    except ValueError:
      try:
        cols[i] = sdic[cols[i].upper()]
      except KeyError:
        lerror(msg + ' (check the roman numeral or use ' + \
               'an integer if stage above 25)',ofolder,verbosity)
        return False, []
      except:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Lower level
    i = 2
    msg = 'Lower level energy '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])*1e-5
    except:
      lerror(msg,ofolder,verbosity)
      return False, []
    i = 3
    msg = 'Lower level angular momentum '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []
    i = 4
    msg = 'Lower level Lande factor '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Upper level
    i = 5
    msg = 'Upper level energy '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])*1e-5
    except:
      lerror(msg,ofolder,verbosity)
      return False, []
    i = 6
    msg = 'Upper level angular momentum '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []
    i = 7
    msg = 'Upper level Lande factor '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Einstein coefficient
    i = 8
    msg = 'Einstein coefficient '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])*1e-8
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Transition lines
    bark = ['s','p','d','f']
    bardic = {'s':0,'p':1,'d':2,'f':3}
    broads = ['barklem','unsold','param','gamma']
    broads_dic = {'barklem':'0','unsold':'1','param':'2','gamma':'3'}


    # VdW broadening
    i = 9
    cols[i] = cols[i].lower()
    # Type
    msg = 'Van der Waals broadening type '+cols[i]+ \
          ' is none of the following: '+ \
          'following: {0}, {1}, {2}, or {3}'.format(*broads)
    if cols[i] not in broads:
      lerror(msg,ofolder,verbosity)
      return False, []
    # Check
    if cols[i] == 'barklem':
      changed = False
      cols[i+1] = cols[i+1].lower()
      msg = 'no Barklem data for l {0}'.format(cols[i+1]) + \
            ', changed to Unsold'
      if cols[i+1] not in bark:
        lerror(msg,ofolder,verbosity)
        cols[i+1] = 1e0
        cols[i+2] = 0e0
        cols[i+3] = 1e0
        cols[i+4] = 0e0
        changed = True
      cols[i+3] = cols[i+3].lower()
      msg = 'no Barklem data for l {0}'.format(cols[i+3]) + \
            ', changed to Unsold'
      if cols[i+3] not in bark:
        lerror(msg,ofolder,verbosity)
        cols[i+1] = 1e0
        cols[i+2] = 0e0
        cols[i+3] = 1e0
        cols[i+4] = 0e0
        changed = True
      if not changed:
        cols[i+1] = bardic[cols[i+1]]
        cols[i+3] = bardic[cols[i+3]]
    if cols[i] == 'unsold' or cols[i] == 'param':
      for j in range(1,5):
        msg = 'broadening parameter '+cols[i+j]
        try:
          cols[i+j] = interpret(cols[i+j])
          cols[i+j] = float(cols[i+j])
        except:
            lerror(msg,ofolder,verbosity)
            return False, []
    if cols[i] == 'gamma':
      msg = 'Van der Waals coefficient '+cols[i]
      try:
        cols[i+1] = interpret(cols[i+1])
        cols[i+1] = 10.0**float(cols[i+1])
      except:
        lerror(msg,ofolder,verbosity)
        return False, []
    cols[i] = broads_dic[cols[i]]

    # Stark coefficient
    i = 14
    msg = 'Stark coefficient '+cols[i]
    if '--' in cols[i] or '-+' in cols[i]:
        ff = -1e0
        cols[i] = cols[i][1:]
        if cols[i][0] == '+':
            cols[i] = cols[i][1:]
    else:
        ff = 1e0
    try:
      cols[i] = interpret(cols[i])
      cols[i] = ff*(10.0**(float(cols[i])))
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Collisional coefficient
    i = 15
    msg = 'Collisional oscillator strength '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Radiative broadening
    i = 16
    msg = 'Radiative coefficient '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = 10.0**(float(cols[i]))
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Total frequencies
    i = 17
    msg = 'Number of frequencies '+cols[i]
    try:
      cols[i] = int(cols[i])
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Core frequencies
    i = 18
    msg = 'Number of frequencies for line core '+cols[i]
    try:
      cols[i] = int(cols[i])
      if cols[i] > cols[i-1]:
        msg = 'Number of frequencies for line core '+cols[i]+ \
              'larger than total '+cols[i-1]+','
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Total Doppler widths
    i = 19
    msg = 'Number of Doppler widths '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
      if cols[i] < 0:
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Core Doppler widths
    i = 20
    msg = 'Number of Doppler widths for line core '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
      if cols[i] > cols[i-1]:
        msg = 'Number of Doppler widths for line core '+cols[i]+ \
              'larger than total '+cols[i-1]+','
        lerror(msg,ofolder,verbosity)
        return False, []
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # Line type
    classe = ['e1','m1','e2','m2','un']
    classe_dic = {'e1':1,'m1':2,'e2':3,'m2':4,'un':5}
    i = 21
    msg = 'Type of transition '+cols[i].lower()
    try:
      cols[i] = classe_dic[cols[i].lower()]
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # tau limit
    i = 22
    msg = 'Tau limit '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # T limit
    i = 23
    msg = 'Temperature limit '+cols[i]
    try:
      cols[i] = interpret(cols[i])
      cols[i] = float(cols[i])
    except:
      lerror(msg,ofolder,verbosity)
      return False, []

    # If we are here, everything is valid
    return True, [iatom_type] + cols


  # Routine to process an LTE line input
  def process_LTEline(entries,ofolder,verbosity):
    ''' Process the entries for LTE lines
    '''

    # Initialize
    NL = 0
    out = []

    # For each entry
    for entry in entries:

      # If file
      if os.path.isfile(entry):

          # Open file
          with open(entry) as g:

            # For each line
            for line in g:

              # Manage Comments
              lline = line.strip()
              j = lline.find('!')
              if j != -1:
                lline = lline[:j]
              j = lline.find('#')
              if j != -1:
                lline = lline[:j]
              if len(lline.split()) == 0:
                continue

              # Try as Kurucz
              valid, lout = process_LTEline_entry_kur(lline)

              # Process the entry as custom
              if not valid:
                valid, lout = process_LTEline_entry(lline,ofolder, \
                                                    verbosity)

              # If valid entry, add to output
              if valid:
                NL += 1
                for lou in lout:
                  out.append(lou)

      # Not a file
      else:

        # Try as Kurucz
        valid, lout = process_LTEline_entry_kur(entry)

        # Process the entry
        if not valid:
          valid, lout = process_LTEline_entry(entry,ofolder,verbosity)

        # If valid entry, add to output
        if valid:
          NL += 1
          for lou in lout:
            out.append(lou)

    # Return results
    return NL, out


  #
  # Argument control
  #

  # Requires one argument
  if len(sys.argv) < 1:
    sys.exit(' # At least one argument needed')

  # Try getting ID
  try:
    dni = sys.argv[2]
  except:
    dni = '000000000'

  # Try to open file
  try:
    f=open(sys.argv[1],'r')
  except:
    sys.exit(' # No input file found')

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
    if ((len(line.split()) == 0) or (line.count('=') == 0)):
      continue
    else:
      lines_n.append(line)
  lines = lines_n

  # Start output file
  filename = 'tmp_input_'+dni
  f = open(filename,'w')
  f.write('1\n')

  # Dictionary of fields that are case sensitive
  URL = ['ATOM_INPUT','ATOM_BACK','OUT_FOLDER','ATMO_INPUT', \
         'CONTINUUM_INPUT','BFIELD_INPUT','SOLUTION_INPUT', \
         'MOLECULE_INPUT','OPACITY_FUDGE','KURUCZ', \
         'WAVELENGTHS','ASYMM_INPUT','PARFUN','ABUND','BARK_SP', \
         'BARK_PD','BARK_DF','MPIDETAIL','OPERFORM','CHIANTI_PATH', \
         'SPECT_INPUT','ATOM_ION','ATOM_POPU','ATOM_FIX_POP', \
         'ATOM_ZERO_ION','DATA_FILE','INV_INIT','ATOM_NO_WAVE', \
         'ATOM_FIX_POP_LTERM','LTE_LINE','WEIGHT_FILE','PSF_FWHM', \
         'INV_MASK']

  # Dictionary of fields that are additive
  APP = ['ATOM_INPUT','ATOM_BACK','MOLECULE_INPUT','KURUCZ', \
         'ASYMM_INPUT','ATOM_FIX_POP','ATOM_ZERO_ION', \
         'ATOM_POPU', 'ATOM_ION','ATOM_FIX_POP_LTERM', \
         'LIM_STK','LIM_CTR','LIM_TAU', \
         'LIM_COLS_TT','LIM_COLS_LL','LIM_DAMP','LIM_BACK', \
         'LIM_POP','LIM_QEL','ATMO_STRAT','WEIGHT','ATOM_NO_WAVE', \
         'PSF_FWHM','LTE_LINE','K_CUT_TERM','EXCLUDE_PIXEL', \
         'WEIGHT_FACTOR','SIGMA_FACTOR']

  # Inversion variables
  varis = ['B','BT','BP','F','T','VX','VY','VZ','VT','PG', \
           'J21R','J21I','J22R','J22I']

  # Hidden inversion variables
  varis_hidden = ['J21R','J21I','J22R','J22I']

  # Add to APP the especial boundary keywords
  for var in varis:

      # Except for the diffuse light factor
      if var == 'F': continue

      # Append for inversion variables
      APP.append('EBOUNDS_'+var)

  # Initialize reading dictionary
  Dictionary = {}

  # For each line
  for line in lines:

    # Get label and value
    Key, Target = line.split('=', 1)

    # Check if URL (case sensitive)
    url_bool = False
    for url in URL:
      if Key.strip().upper() == url:
        url_bool = True
        break

    # Check if APP (additive)
    app_bool = False
    for app in APP:
      if Key.strip().upper() == app:
        app_bool = True
        break

    # If URL
    if url_bool:

      # If APP
      if app_bool:

        # Append if already exists
        if Key.strip().upper() in Dictionary:
          Dictionary[Key.strip().upper()].append(Target.strip())
        # Initialize otherwise
        else:
          Dictionary[Key.strip().upper()] = [Target.strip()]

      # If not APP
      else:

        # Just write
        Dictionary[Key.strip().upper()] = (Target.strip()).split()

    # Not an URL
    else:

      # If APP
      if app_bool:

        # Append if already exists
        if Key.strip().upper() in Dictionary:
          Dictionary[Key.strip().upper()].append( \
                                   (Target.strip().upper()).split())
        # Initialize otherwise
        else:
          Dictionary[Key.strip().upper()] = \
                                  [(Target.strip().upper()).split()]

      # If not APP
      else:

        # Just write
        Dictionary[Key.strip().upper()] = \
                                   (Target.strip().upper()).split()

  #
  # Make uppercase the non URL
  #

  # For each key
  for key in list(Dictionary.keys()):

    # Initialize flags
    url_bool = False
    app_bool = False

    # For each case-sensitive key
    for url in URL:

      # If same key
      if key == url:

        # Found case-sensitive
        url_bool = True
        break

    # For each additive key
    for app in APP:

      # If same key
      if key == app:

        # Found additive
        app_bool = True
        break

    # If no case-sensitive
    if not url_bool:

      # If additive
      if app_bool:

        # For each value in dictionary
        for vals in Dictionary[key]:

          # For each value, make uppercase
          for val in vals: val = val.upper()

      # Not additive
      else:

        # For each value in dictionary, make uppercase
        for val in Dictionary[key]: val = val.upper()

  ###################################################################
  # HANLERT
  ###################################################################

  #
  # Preliminar checks
  # Check verbosity
  #

  # VERBOSE
  check = 0
  if 'VERBOSE' in Dictionary:
    val = Dictionary['VERBOSE'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      check = 1
      verbosity = True
    if val == 'N' or val == 'NO' or val == 'NON':
      check = 1
      verbosity = False
  if check == 0:
      verbosity = False

  # OUTFILE
  if 'OUT_FOLDER' in Dictionary:
    val = Dictionary['OUT_FOLDER'][0]
    ofolder = val
  else:
    ofolder = 'Outputs/Default'

  # RUN_MODE
  check = 0
  if 'RUN_MODE' in Dictionary:
    val = Dictionary['RUN_MODE'][0]
    if '1D' in val and 'S' in val:
        rmode = 0
    elif ('1.5D' in val or '15D' in val) and 'S' in val:
        rmode = 1
    elif 'CLE' in val:
        rmode = 2
    elif 'INV' in val:
        rmode = -1
    else:
        rmode = 0
  else:
    rmode = 0

  # SKIP_FILE
  # To know if we need solution folder
  if 'SKIP_SFILE' in Dictionary:
    val = Dictionary['SKIP_SFILE'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      no_folder = True
    else:
      no_folder = False
  else:
    no_folder = False
  if rmode >= 0:
    if 'KEEP_RHOKQ' in Dictionary:
      val = Dictionary['KEEP_RHOKQ'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        no_folder = False
  if rmode >= 0:
    if 'KEEP_JKQ' in Dictionary:
      val = Dictionary['KEEP_JKQ'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        no_folder = False
  if rmode >= 0:
    if 'KEEP_STOKES_QUAD' in Dictionary:
      val = Dictionary['KEEP_STOKES_QUAD'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        no_folder = False
  if rmode >= 0:
    if 'KEEP_JKQNU' in Dictionary:
      val = Dictionary['KEEP_JKQNU'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        no_folder = False

  # Check the folder and try to create it if does not exits
  exist = os.path.isdir(ofolder)
  if not exist:
    try:
      os.makedirs(ofolder)
      verbose(' - Output folder did not exist, I created ' + \
              'it for you', ofolder, verbosity)
    except:
      verbose(' # Could not create output folder', '', verbosity)
      abort(f, filename)
  elif not(verbosity):
    exist = os.path.isfile(ofolder+'/verbose')
    if exist:
      os.remove(ofolder+'/verbose')

  # FILE BUILD

  # RUN_MODE
  check = 0
  if 'RUN_MODE' in Dictionary:
    val = Dictionary['RUN_MODE'][0]
    if '1D' in val and 'S' in val:
        f.write('0\n')
        rmode = 0
    elif ('1.5D' in val or '15D' in val) and 'S' in val:
        f.write('1\n')
        rmode = 1
    elif 'CLE' in val:
        f.write('2\n')
        rmode = 2
    elif 'INV' in val:
        f.write('-1\n')
        rmode = -1
    else:
        verbose(' # RUN_MODE not recognized keyword', \
                ofolder, verbosity)
        abort(f, filename)
  else:
    verbose(' # RUN_MODE necessary keyword', ofolder, verbosity)
    abort(f, filename)

  # ATMO_INPUT
  if rmode >= 0:
    if 'ATMO_INPUT' in Dictionary:
      val = Dictionary['ATMO_INPUT'][0]
      f.write(val+'\n')
      try:
        f2=open(val)
        f2.close()
      except:
        verbose(' # ATMO_INPUT file not found '+val, \
                ofolder, verbosity)
        abort(f, filename)
    else:
      verbose(' # ATMO_INPUT necessary keyword', ofolder, verbosity)
      abort(f, filename)
  elif rmode == -1:
    f.write('NONE\n')

  # DATA_FILE
  if rmode < 0:
    if 'DATA_FILE' in Dictionary:
      val = Dictionary['DATA_FILE'][0]
      try:
        f2=open(val)
        f2.close()
      except:
        verbose(' # DATA_FILE file not found '+val, ofolder, verbosity)
        abort(f, filename)
    else:
      verbose(' # DATA_FILE necessary keyword', ofolder, verbosity)
      abort(f, filename)

  # ATMO_SCALE
  if rmode == 1 or rmode == -1:
    if 'ATMO_SCALE' in Dictionary:
      val = Dictionary['ATMO_SCALE'][0]
      scale = val.split()
      if len(scale) >= 2:
        try:
          scale[1] = interpret(scale[1])
          tf = '{0}\n'.format(1e2/float(scale[1]))
        except ValueError:
          verbose(' # ATMO_SCALE wavelength seems to not be float', \
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
      else:
        verbose(' # ATMO_SCALE scale not recognized', \
                verbfile, verbosity)
        abort(f,filename)
      f.write(tf)
    else:
      f.write('H\n')
      f.write('0.2\n')
  else:
    f.write('N\n')
    f.write('0.2\n')

  # ATMO_CHAR
  if rmode == 1 or rmode == 2 or rmode == -1:
    if 'ATMO_CHAR' in Dictionary:
      val = Dictionary['ATMO_CHAR'][0]
      # If specifying number densities
      if 'NE' in val and 'NH' in val:
        f.write('0\n')
      # If specifying only electron number density
      elif 'NE' in val:
        f.write('1\n')
      # If specifying only electron pressure
      elif 'PE' in val:
          f.write('2\n')
      # If specifying only electron mass density
      elif 'RHOE' in val:
          if rmode != 0:
              verbose(' # ATMO_CHAR option "rhoe"' + \
                      ' has been deprecated',verbfile,verbosity)
              abort(f,filename)
          f.write('3\n')
      # If specifying only gas pressure
      elif 'PG' in val:
          f.write('4\n')
      # If specifying only gas density
      elif 'RHO' in val:
          f.write('5\n')
      # If we expect number densities
      else:
          f.write('0\n')
    # Default number densities
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # Initialize atoms output buffer, labels, and mapping
  atom_lab_act = []
  atom_lab_pas = []
  atom_map = {}

  # ATOM_INPUT
  if 'ATOM_INPUT' in Dictionary:

    # Get atom entries
    val = Dictionary['ATOM_INPUT']

    # Initialize number of active atoms and add to buffer
    NA = len(val)

    # For each name
    for iname,name in enumerate(val):

      # Split in elements
      cols = name.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
          if cols[ii] == '': cols.pop(ii)

      # If only one column, label is number
      if len(cols) < 2:
        label = '{0}'.format(iname)
      # If two columns, label is second column
      else:
        label = '{0}'.format(cols[1].strip())

      # Check repeated
      if label in atom_lab_act:
        verbose(' # ATOM_INPUT label repeated'+label, \
                ofolder, verbosity)
        abort(f, filename)

      # Add label to list
      atom_lab_act.append(label)

      # Initialize atom entry
      atom_map[label] = [cols[0]+'\n','N\n','N\n','N\n', \
                         '-1\n','N\n','N\n']

      # Check if atomic model file exists
      try:
        f2=open(cols[0])
        f2.close()
      except:
        verbose(' # ATOM_INPUT file not found '+cols[0], \
                ofolder, verbosity)
        abort(f, filename)
  else:
    verbose(' # ATOM_INPUT necessary keyword', ofolder, verbosity)
    abort(f, filename)


  # ATOM_BACK
  if 'ATOM_BACK' in Dictionary:

    # Get background atom entries
    val = Dictionary['ATOM_BACK']

    # Initialize number of passive atoms and add to buffer
    NAb = len(val)

    # For each name
    for jname,name in enumerate(val):

      # True label
      iname = jname + NA

      # Split in elements
      cols = name.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
        if cols[ii] == '': cols.pop(ii)

      # If only one column, label is number
      if (len(cols) < 2):
        label = '{0}'.format(iname)
      # If two columns, label is second
      else:
        label = '{0}'.format(cols[1].strip())

      # Check repeated
      if label in atom_lab_act or label in atom_lab_pas:
        verbose(' # ATOM_BACK label repeated'+label, \
                ofolder, verbosity)
        abort(f, filename)

      # Add label to list
      atom_lab_pas.append(label)

      # Initialize atom entry
      atom_map[label] = [cols[0]+'\n','N\n']

      # Check if atomic model exists
      try:
        f2=open(cols[0])
        f2.close()
      except:
        verbose(' # ATOM_BACK file not found '+atom, \
                ofolder, verbosity)
        abort(f, filename)
  else:
    verbose(' # Set ATOM_BACK=NONE(DEFAULT)', ofolder, verbosity)
    NAb = 0

  # ATOM_POPU
  if 'ATOM_POPU' in Dictionary:

    # Get populations entries
    val = Dictionary['ATOM_POPU']

    # For each name
    for iname,name in enumerate(val):

      # Split in elements
      cols = name.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
        if cols[ii] == '': cols.pop(ii)

      # If only one column, label is number
      if len(cols) < 2:
        label = '{0}'.format(iname)
      # If two columns, label is second
      else:
        label = cols[1].strip()

      # If label does not exists, problem
      if label not in atom_lab_act and \
         label not in atom_lab_pas:
        verbose(' # ATOM_POPU Label not found among atoms '+label, \
                ofolder, verbosity)
        abort(f, filename)

      # Add to atom entry
      atom_map[label][1] = cols[0]+'\n'

      # Check if atomic model file exists
      try:
        f2=open(cols[0])
        f2.close()
      except:
        verbose(' # ATOM_POPU file not found '+atom, \
                ofolder, verbosity)
        abort(f, filename)

  # ATOM_FIX_POP
  if 'ATOM_FIX_POP' in Dictionary:

    # Get fix populations entries
    val = Dictionary['ATOM_FIX_POP']

    # For each label
    for labels in val:

      # Try to split
      cols = labels.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
        if cols[ii] == '': cols.pop(ii)

      # For each label in each entry
      for label in cols:

        # If label does not exists, problem
        if label not in atom_lab_act:
          verbose(' # ATOM_FIX_POP Label not found ' + \
                  'among atoms '+label, ofolder, verbosity)
          abort(f, filename)

        # Set atom to fix population
        atom_map[label][2] = 'F\n'

  # ATOM_ZERO_ION
  if 'ATOM_ZERO_ION' in Dictionary:

    # Get fix populations entries
    val = Dictionary['ATOM_ZERO_ION']

    # For each label
    for labels in val:

      # Try to split
      cols = labels.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
        if cols[ii] == '': cols.pop(ii)

      # For each label in each entry
      for label in cols:

        # If label does not exists, problem
        if label not in atom_lab_act:
          verbose(' # ATOM_ZERO_ION Label not found ' + \
                  'among atoms '+label, ofolder, verbosity)
          abort(f, filename)

        # Set atom to nullify last ion
        atom_map[label][3] = 'F\n'

  # ATOM_NO_WAVE
  if 'ATOM_NO_WAVE' in Dictionary:

    # Get fix populations entries
    val = Dictionary['ATOM_NO_WAVE']

    # For each label
    for labels in val:

      # Try to split
      cols = labels.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
        if cols[ii] == '': cols.pop(ii)

      # For each label in each entry
      for label in cols:

        # If label does not exists, problem
        if label not in atom_lab_act:
          verbose(' # ATOM_NO_WAVE Label not found ' + \
                  'among atoms '+label, ofolder, verbosity)
          abort(f, filename)

        # Set atom to fix population
        atom_map[label][5] = 'Y\n'

  # ATOM_FIX_POP_LTERM
  if 'ATOM_FIX_POP_LTERM' in Dictionary:

    # Get fix populations entries
    val = Dictionary['ATOM_FIX_POP_LTERM']

    # For each label
    for labels in val:

      # Try to split
      cols = labels.split(' ')

      # Remove additional spaces
      for ii in range(len(cols)-1,-1,-1):
        if cols[ii] == '': cols.pop(ii)

      # For each label in each entry
      for label in cols:

        # If label does not exists, problem
        if label not in atom_lab_act:
          verbose(' # ATOM_FIX_POP_LTERM Label not found ' + \
                  'among atoms '+label, ofolder, verbosity)
          abort(f, filename)

        # Set atom to fix population
        atom_map[label][6] = 'F\n'

  # ATOM_ION
  if rmode == 2:

    # Only for CLE
    if 'ATOM_ION' in Dictionary:

      # Get values
      val = Dictionary['ATOM_ION']

      # For each entry
      for iname,name in enumerate(val):

        # Split in elements
        cols = name.split(' ')

        # Remove additional spaces
        for ii in range(len(cols)-1,-1,-1):
          if cols[ii] == '': cols.pop(ii)

        # If only one column, label is number
        if len(cols) < 2:
          label = '{0}'.format(iname)
        # If two columns, label is second
        else:
          label = cols[1].strip()

        # If label does not exists, problem
        if label not in atom_lab_act and \
           label not in atom_lab_pas:
          verbose(' # ATOM_ION Label not found among atoms '+label, \
                  ofolder, verbosity)
          abort(f, filename)

        # If float
        check = 0
        try:
          cols[0] = interpret(cols[0])
          num = float(cols[0])
          # Update atom entry
          atom_map[label][4] = '1\n{0}\n'.format(num)
          check = 1
        except:
          pass

        # Not float, could be path
        if check == 0:
          try:
            f2=open(cols[0])
            f2.close()
            # Update atom entry
            atom_map[label][4] = '0\n'+cols[0]+'\n'
          except:
            msg = ' # ATOM_ION Problem reading ionization file {0}'
            verbose(msg.format(cols[0]),ofolder, verbosity)
            abort(f,filename)

  #
  # Now we can write the atomic data
  #

  # Number of active atoms
  f.write('{0}\n'.format(NA))

  # For each active atom
  for label in atom_lab_act:
    # Write the data in the entry
    for data in atom_map[label]:
      f.write(data)

  # Number of passive atoms
  f.write('{0}\n'.format(NAb))

  # For each active atom
  for label in atom_lab_pas:
    # Write the data in the entry
    for data in atom_map[label]:
      f.write(data)

  # MOLECULE_INPUT
  if 'MOLECULE_INPUT' in Dictionary:
    val = Dictionary['MOLECULE_INPUT']
    NM = len(val)
    f.write(str(NM)+'\n')
    for name in val:
      f.write(name+'\n')
      try:
        f2=open(name)
        f2.close()
      except:
        verbose(' # MOLECULE_INPUT file not found '+name, \
                ofolder, verbosity)
        abort(f, filename)
  else:
    verbose(' # Set MOLECULE_INPUT=NONE(DEFAULT)', ofolder, verbosity)
    f.write('0\n')

  # BFIELD_INPUT
  if rmode == 0 or rmode == -1:
    check = 0
    if 'BFIELD_INPUT' in Dictionary:
      val = Dictionary['BFIELD_INPUT']
      ncol = len(val)
      # If three elements
      if ncol == 3:
        bb = []
        for b in val:
          try:
            b = interpret(b)
            bb.append(float(b))
          except:
            break
        if len(bb) == 3:
          f.write('N\n')
          for b in val:
            f.write(b+'\n')
          check = 1
      else:
        val = val[0]
        try:
          f2=open(val)
          f2.close()
          f.write('F\n')
          f.write(val+'\n')
          check = 1
        except:
          f.write('F\n')
          f.write('NONE\n')
          check = 1
    if check == 0:
      f.write('F\n')
      f.write('NONE\n')
  else:
    f.write('F\n')
    f.write('NONE\n')

  # OPACITY_FUDGE
  if 'OPACITY_FUDGE' in Dictionary:
    val = Dictionary['OPACITY_FUDGE'][0]
    try:
      f2=open(val)
      f2.close()
      f.write(val+'\n')
    except:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # SPECT_INPUT
  if rmode == 2:
    if 'SPECT_INPUT' in Dictionary:
      val = Dictionary['SPECT_INPUT'][0]
      try:
        f2=open(val)
        f2.close()
        f.write(val+'\n')
      except:
        f.write('NONE\n')
    else:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # USE_ALLEN
  if rmode == 2:
    check = 0
    if 'USE_ALLEN' in Dictionary:
      val = Dictionary['USE_ALLEN'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # FLAT_CLE_IN
  if rmode == 2:
    check = 0
    if 'FLAT_CLE_IN' in Dictionary:
      val = Dictionary['FLAT_CLE_IN'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # CHIANTI_PATH
  if rmode == 2:
    if 'CHIANTI_PATH' in Dictionary:
      val = Dictionary['CHIANTI_PATH'][0]
      try:
        if os.path.isdir(val):
          if val[-1] != '/': val += '/'
          f.write(val+'\n')
        else:
          f.write('NONE\n')
      except:
        f.write('NONE\n')
    else:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # T_RAD
  if rmode == 2:
    check = 0
    if 'T_RAD' in Dictionary:
      val = Dictionary['T_RAD'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except ValueError:
        pass
      except:
        raise
    if check == 0:
      f.write('5d3\n')
  else:
    f.write('0d0\n')

  # R_STAR
  if rmode == 2:
    check = 0
    if 'R_STAR' in Dictionary:
      val = Dictionary['R_STAR'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except ValueError:
        pass
      except:
        raise
    if check == 0:
      f.write('6.95700d10\n')
  else:
    f.write('0d0\n')

  # NEGLECT_CONTINUUM
  if rmode == 2:
    check = 0
    if 'NEGLECT_CONTINUUM' in Dictionary:
      val = Dictionary['NEGLECT_CONTINUUM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # PARFUN
  if 'PARFUN' in Dictionary:
    val = Dictionary['PARFUN'][0]
    try:
      f2=open(val)
      f2.close()
      f.write(val+'\n')
    except:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # ABUND
  if 'ABUND' in Dictionary:
    val = Dictionary['ABUND'][0]
    try:
      f2=open(val)
      f2.close()
      f.write(val+'\n')
    except:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # BARK_SP
  if 'BARK_SP' in Dictionary:
    val = Dictionary['BARK_SP'][0]
    try:
      f2=open(val)
      f2.close()
      f.write(val+'\n')
    except:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # BARK_PD
  if 'BARK_PD' in Dictionary:
    val = Dictionary['BARK_PD'][0]
    try:
      f2=open(val)
      f2.close()
      f.write(val+'\n')
    except:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # BARK_DF
  if 'BARK_DF' in Dictionary:
    val = Dictionary['BARK_DF'][0]
    try:
      f2=open(val)
      f2.close()
      f.write(val+'\n')
    except:
      f.write('NONE\n')
  else:
    f.write('NONE\n')

  # IGNORE_BB
  check = 0
  # if CLE, always ignore
  if rmode == 2:
    f.write('Y\n')
  else:
    if 'IGNORE_BB' in Dictionary:
      val = Dictionary['IGNORE_BB'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

  # KURUCZ
  if 'KURUCZ' in Dictionary:
    val = Dictionary['KURUCZ']
    NK = 0
    out = []
    for name in val:
      try:
        f2=open(name)
        f2.close()
        NK += 1
        out.append(name)
      except:
        pass
    f.write(str(NK)+'\n')
    for ou in out:
      f.write(ou+'\n')
  else:
    f.write('0\n')

  # LTE_LINE
  if rmode != 2:
    if 'LTE_LINE' in Dictionary:
      vals = Dictionary['LTE_LINE']
      NL, out = process_LTEline(vals,ofolder,verbosity)
      f.write(str(NL)+'\n')
      for ou in out:
        f.write(str(ou)+'\n')
    else:
      f.write('0\n')
  # No for CLE
  else:
    f.write('0\n')

  # WAVELENTHS
  if 'WAVELENGTHS' in Dictionary:
    val = Dictionary['WAVELENGTHS']
    NK = 0
    out = []
    for name in val:
      try:
        f2=open(name)
        f2.close()
        NK += 1
        out.append(name)
      except:
        pass
    f.write(str(NK)+'\n')
    for ou in out:
      f.write(ou+'\n')
  else:
    f.write('0\n')

  # ASYMM_INPUT
  if rmode == 0:
    if 'ASYMM_INPUT' in Dictionary:
      val = Dictionary['ASYMM_INPUT']
      out = []
      nentry = 0
      nentryV = 0
      nentryF = 0
      passed = True
      for lval in val:
        if lval.strip().lower() == 'none':
          continue
        try:
          cols = lval.split()
          if len(cols) == 4:
            K = int(cols[0])
            Q = int(cols[1])
            cols[2] = interpret(cols[2])
            cols[3] = interpret(cols[3])
            ar = float(cols[2])
            ai = float(cols[3])
            if K < 1 or K > 2:
              passed = False
              verbose(' # Input asymmetry must have ' + \
                      'multipole K in the [1,2] range, '+ \
                      'your input is {0}'.format(K), \
                      ofolder,verbosity)
              break
            if Q < 0 or Q > K:
              passed = False
              verbose(' # Input asymmetry must have ' + \
                      'multipole Q in the [0,K] range, ' + \
                      'your input is {0}'.format(Q), \
                      ofolder,verbosity)
              break
            if (ar*ar + ai*ai) <= 0. or math.sqrt(ar*ar + ai*ai) > 1.:
              passed = False
              verbose(' # Input asymmetry absolute value must ' + \
                      'be in the (0,1] range, your input is ' + \
                      '{0} + i{1}'.format(ar,ai), ofolder, verbosity)
              break
            nentry += 1
            nentryV += 1
            out.append('V')
            out.append('{0} {1} {2} {3}'.format(K,Q,ar,ai))
          elif len(cols) == 1:
            fil = lval
            try:
              f2=open(fil)
              f2.close()
              out.append('F')
              out.append(fil)
              nentry += 1
              nentryF += 1
            except IOError:
              passed = False
              verbose(' # ASYMM_INPUT file not found '+fil, \
                      ofolder, verbosity)
              break
            except:
              raise
          else:
            passed = False
            verbose(' # Number of elements in ASYMM_INPUT ' + \
                    'not equal to 1 or 4 in line: ' + lval + ',', \
                    ofolder, verbosity)
            break
        except:
          raise
      if nentry > 0:
        if passed:
          out = ['{0}'.format(nentry),'{0}'.format(nentryV), \
                 '{0}'.format(nentryF)] + out
          for o in out:
            f.write(o+'\n')
      else:
        f.write('0\n')
    else:
      f.write('0\n')
  # 1.5D admits file
  elif rmode == 1 or rmode == -1:
    if 'ASYMM_INPUT' in Dictionary:
      val = Dictionary['ASYMM_INPUT'][-1]
      if val.strip().lower() != 'none':
        fil = val.strip()
        try:
          f2=open(fil)
          f2.close()
          f.write('1\n0\n1\nF\n{0}\n'.format(fil))
          passed = True
        except IOError:
          passed = False
          verbose(' # ASYMM_INPUT file not found '+fil, \
                  ofolder, verbosity)
        except:
          raise
        if not passed:
          f.write('0\n')
      else:
        f.write('0\n')
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # FORCE_ASYMM
  check = 0
  if 'FORCE_ASYMM' in Dictionary:
    val = Dictionary['FORCE_ASYMM'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # STIM
  check = 0
  if 'STIM' in Dictionary:
    val = Dictionary['STIM'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('Y\n')

  # OUTFILE
  if 'OUT_FOLDER' in Dictionary:
    val = Dictionary['OUT_FOLDER'][0]
    f.write(val+'\n')
  else:
    f.write('Outputs/Default\n')

  # POLAR_NODES
  if rmode >= -1 and rmode <= 2:
    check = 0
    if 'POLAR_NODES' in Dictionary:
      val = Dictionary['POLAR_NODES'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
        pnodes = int(val)
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(4))
      pnodes = 4
  else:
    f.write('0\n')
    pnodes = 0

  # AXIAL_NODES
  if rmode >= -1 and rmode <= 2:
    check = 0
    if 'AXIAL_NODES' in Dictionary:
      val = Dictionary['AXIAL_NODES'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
        anodes = int(val)
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(2))
      anodes = 2
  else:
    f.write('0\n')
    anodes = 0

  # POLARI_NODES
  if rmode >= -1 and rmode <= 2:
    check = 0
    if 'POLARI_NODES' in Dictionary:
      val = Dictionary['POLARI_NODES'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(pnodes))
  else:
    f.write('0\n')

  # AXIALI_NODES
  if rmode >= -1 and rmode <= 2:
    check = 0
    if 'AXIALI_NODES' in Dictionary:
      val = Dictionary['AXIALI_NODES'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(anodes))
  else:
    f.write('0\n')

  # POLAR_LOS
  if rmode == 0 or rmode == 1:
    check = 0
    if 'POLAR_LOS' in Dictionary:
      val = Dictionary['POLAR_LOS']
      out = []
      for k in range(len(val)):
        val[k] = interpret(val[k])
        try:
          out.append(float(val[k]))
        except ValueError:
          try:
            verbose(' # No float in mu LOS: ',val[k], ofolder, \
                    verbosity)
          except:
            pass
        except:
          try:
            verbose(' # Failed writting mu LOS: ',val[k], ofolder, \
                    verbosity)
          except:
            pass
      if len(out) > 0:
        f.write('{0:7d}\n'.format(len(out)))
        for ou in out:
          f.write('{0:22.16e} '.format(ou))
        f.write('\n')
        check = 1
        nopol = False
    if check == 0:
      f.write('0\n')
      nopol = True
  elif rmode == 2:
    f.write('1\n0.\n')
  elif rmode == -1:
    f.write('1\n0.\n')
  else:
    f.write('0\n')

  # AXIAL_LOS
  if rmode == 0 or rmode == 1:
    check = 0
    if 'AXIAL_LOS' in Dictionary:
      val = Dictionary['AXIAL_LOS']
      out = []
      if nopol:
        verbose(' # Cannot specify AXIAL_LOS without POLAR_LOS', \
                ofolder, verbosity)
        abort(f, filename)
      for k in range(len(val)):
        val[k] = interpret(val[k])
        try:
          out.append(float(val[k]))
        except ValueError:
          try:
            verbose(' # No float in azimuth LOS: ',val[k], ofolder, \
                    verbosity)
          except:
            pass
        except:
          try:
              verbose(' # Failed writting azimuth LOS: ',val[k], \
                      ofolder, verbosity)
          except:
            pass
      if len(out) > 0:
        f.write('{0:7d}\n'.format(len(out)))
        for ou in out:
          f.write('{0:22.16e} '.format(ou))
        f.write('\n')
        check = 1
    if check == 0:
      if nopol:
        f.write('0\n')
      else:
        f.write('{0:7d}\n{1:22.16e} \n'.format(1,0))
  elif rmode == 2:
    f.write('1\n0.\n')
  elif rmode == -1:
    f.write('1\n0.\n')
  else:
    f.write('0\n')

  # MODE
  if rmode >= 0:
    check = 0
    if 'MODE' in Dictionary:
      val = Dictionary['MODE'][0]
      if val == 'R' or val == 'RE' or val == 'REA' or val == 'READ':
        f.write('R\n')
        MOD = 'R'
        check = 1
      if val == 'S' or val == 'SO' or val == 'SOL' or val == 'SOLV' \
         or val == 'SOLVE':
        f.write('W\n')
        MOD = 'W'
        check = 1
      if val == 'B' or val == 'BO' or val == 'BOT' or val == 'BOTH':
        f.write('B\n')
        MOD = 'B'
        check = 1
    if check == 0:
      f.write('W\n')
      MOD = 'W'
  else:
    f.write('W\n')
    MOD = 'W'

  # FORCE
  if rmode >= 0:
    check = 0
    if 'FORCE' in Dictionary:
      val = Dictionary['FORCE'][0]
      if val == 'I' or val == 'IN' or val == 'INT' or \
         val == 'INTE' or val == 'INTEN' or val == 'INTENS' or \
         val == 'INTENSI' or val == 'INTENSIT' or \
         val == 'INTENSITY' or val == 'OI':
        f.write('I\n')
        check = 1
      if val == 'P' or val == 'PO' or val == 'POL' or \
         val == 'POLA' or val == 'POLAR' or val == 'POLARI' or \
         val == 'POLARIZ' or val == 'POLARIZA' or \
         val == 'POLARIZAT' or val == 'POLARIZATI' or \
         val == 'POLARIZATIO' or val == 'POLARIZATION' or val == 'OP':
        f.write('P\n')
        check = 1
      if val == 'A' or val == 'AL' or val == 'ALL':
        f.write('A\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON' or val == 'NONE':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RESTRICT_TAUC_STRICT
  if rmode == 0 or rmode == 1:
    check = 0
    if 'RESTRICT_TAUC_STRICT' in Dictionary:
      val = Dictionary['RESTRICT_TAUC_STRICT'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RESTRICT_TAUC
  if rmode == 0 or rmode == 1:
    if 'RESTRICT_TAUC' in Dictionary:
      val = Dictionary['RESTRICT_TAUC']
      try:
        val[0] = interpret(val[0])
        val[1] = interpret(val[1])
        t0 = float(val[0])
        t1 = float(val[1])
        t0 = 10e0**t0
        t1 = 10e0**t1
        if t1 < t0:
          ia = t1
          t1 = t0
          t0 = ia
      except:
        verbose(' # RESTRICT_TAUC wrong format', ofolder, verbosity)
        abort(f, filename)
      f.write('Y\n')
      f.write('{0}\n'.format(t0))
      f.write('{0}\n'.format(t1))
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # RESTRICT_HEIGHT_STRICT
  if rmode == 0 or rmode == 1:
    if 'RESTRICT_HEIGHT_STRICT' in Dictionary:
      val = Dictionary['RESTRICT_HEIGHT_STRICT'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RESTRICT_HEIGHT
  if rmode == 0 or rmode == 1:
    if 'RESTRICT_HEIGHT' in Dictionary:
      val = Dictionary['RESTRICT_HEIGHT']
      try:
        val[0] = interpret(val[0])
        val[1] = interpret(val[1])
        t0 = float(val[0])
        t1 = float(val[1])
        if t1 < t0:
          ia = t1
          t1 = t0
          t0 = ia
      except:
        verbose(' # RESTRICT_HEIGHT wrong format', ofolder, verbosity)
        abort(f, filename)
      f.write('Y\n')
      f.write('{0}\n'.format(t0))
      f.write('{0}\n'.format(t1))
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # ZEEMAN_MODE
  check = 0
  zeetypes = {'FULL':'0', 'NOZEEMAN': '1', \
              'NOSPLIT': '2', 'LINEAR': '3'}
  if 'ZEEMAN_MODE' in Dictionary:
    val = Dictionary['ZEEMAN_MODE'][0]
    for ty in list(zeetypes.keys()):
        if val.upper() in ty:
            f.write(zeetypes[ty]+'\n')
            check = 1
            break
  if check == 0:
    f.write('0\n')

  # P_CORR
  check = 0
  if 'P_CORR' in Dictionary:
    val = Dictionary['P_CORR'][0]
    if val == 'N' or val == 'NO':
      f.write('N\n')
      check = 1
    else:
      f.write('Y\n')
      check = 1
  if check == 0:
    f.write('Y\n')

  # FCOL_TRANSFER
  check = 0
  if 'FCOL_TRANSFER' in Dictionary:
    val = Dictionary['FCOL_TRANSFER'][0]
    if val == 'N' or val == 'NO':
      f.write('N\n')
      check = 1
    else:
      f.write('Y\n')
      check = 1
  if check == 0:
    f.write('Y\n')

  # MIT_OFF
  check = 0
  if 'MIT_OFF' in Dictionary:
    val = Dictionary['MIT_OFF'][0]
    if val == 'Y' or val == 'YE' or val == 'YES':
      f.write('-1\n')
      check = 1
    else:
      f.write('1\n')
      check = 1
  if check == 0:
    f.write('0\n')

  # MIT_NODE
  check = 0
  if 'MIT_NODE' in Dictionary:
    val = Dictionary['MIT_NODE'][0]
    try:
      val = interpret(val)
      num = float(val)
      f.write(val+'\n')
      check = 1
    except:
      f.write('1.0\n')
      check = 1
  if check == 0:
    f.write('1.0\n')

  # SKIP_SFILE
  if 'SKIP_SFILE' in Dictionary:
    val = Dictionary['SKIP_SFILE'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # SOLUTION_INPUT
  if rmode == 0:
    if 'SOLUTION_INPUT' in Dictionary:
      val = Dictionary['SOLUTION_INPUT'][0]
      try:
        f2=open(val)
        f2.close()
        f.write(val+'\n')
      except:
        val = ofolder+'/Solution'
        if MOD == 'R' or MOD == 'B':
          try:
            f2=open(val)
            f2.close()
            verbose(' # Set SOLUTION_INPUT='+ofolder+ \
                    '/Solution(DEFAULT)', ofolder, verbosity)
          except:
            f.close()
            verbose(' # Solution file not found', ofolder, \
                    verbosity)
            abort(f, filename)
        f.write(val+'\n')
    else:
      val = ofolder+'/Solution'
      if MOD == 'R' or MOD == 'B':
        try:
          f2=open(val)
          f2.close()
        except:
          f.close()
          verbose(' # Solution file not found', ofolder, verbosity)
          abort(f, filename)
      f.write(val+'\n')
  elif rmode == 1:
    if 'SOLUTION_INPUT' in Dictionary:
      val = Dictionary['SOLUTION_INPUT'][0]
      if list(val)[-1] == '/':
        val = ''.join(list(val)[:-1])
      if os.path.isdir(val):
        f.write(val+'\n')
      else:
        val = ofolder+'/Solution-folder'
        if MOD == 'R' or MOD == 'B':
          if not os.path.isdir(val):
            verbose(' # Solution folder not found', ofolder, \
                    verbosity)
            abort(f, filename)
        f.write(val+'\n')
    else:
      val = ofolder+'/Solution-folder'
      if MOD == 'R' or MOD == 'B':
        if not os.path.isdir(val):
          verbose(' # Solution folder not found', ofolder, \
                  verbosity)
          abort(f, filename)
      f.write(val+'\n')
    # Check the folder and try to create it if does not exits
    exist = os.path.isdir(ofolder+'/Solution-folder')
    if not exist and not no_folder:
      try:
        os.mkdir(ofolder+'/Solution-folder')
        verbose(' - Solution folder did not exist, I created ' + \
                'it for you', ofolder, verbosity)
      except:
        verbose(' # Could not create Solution folder', '', verbosity)
        abort(f, filename)
  elif rmode == 2:
    f.write('NONE\n')
  else:
    val = ofolder+'/Solution'
    f.write(val+'\n')

  # SOLUTION_BACKUP
  if rmode == 0:
    if 'SOLUTION_BACKUP' in Dictionary:
      val = Dictionary['SOLUTION_BACKUP'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        if MOD == 'B' or MOD == 'W':
          found = True
          ind = 0
          name = ofolder+'/Solution'
          name0 = ofolder+'/Solution'
          while (found):
            try:
              f2=open(name)
              f2.close()
              namea = name[::-1]
              ii = namea.find('noituloS')
              name = namea[ii::]
              name = name[::-1]
              name += "_{0:d}".format(ind)
              ind += 1
            except:
              found = False
          if ind > 0:
            try:
              shutil.copyfile(name0,name)
            except:
              verbose(' # Could not backup Solution file', \
                      ofolder, verbosity)
      else:
        f.write('N\n')
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # SOLUTION_KEEPI
  if rmode == 0 or rmode == 1:
    if 'SOLUTION_KEEPI' in Dictionary:
      val = Dictionary['SOLUTION_KEEPI'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
      else:
        f.write('N\n')
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # SOLUTION_KEEPS
  if rmode >= 0:
    if 'SOLUTION_KEEPS' in Dictionary:
      val = Dictionary['SOLUTION_KEEPS'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
      else:
        f.write('N\n')
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # ANISOTROPY_FOCUS
  if rmode < 2:
    if 'ANISOTROPY_FOCUS' in Dictionary:
      val = Dictionary['ANISOTROPY_FOCUS'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
      else:
        f.write('N\n')
    else:
      f.write('N\n')
  else:
    f.write('N\n')

  # K_CUT
  check = 0
  if 'K_CUT' in Dictionary:
    val = Dictionary['K_CUT'][0]
    try:
      f.write('{0:7d}\n'.format(int(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:7d}\n'.format(-1))

  # K_CUTAB
  check = 0
  if 'K_CUTAB' in Dictionary:
    val = Dictionary['K_CUTAB'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # K_CUT_TERM
  check = 0
  if 'K_CUT_TERM' in Dictionary:
    vals = Dictionary['K_CUT_TERM']
    NL = len(vals)
    ranges = []
    for val in vals:
      if len(val) < 3 or len(val) > 4:
        verbose(' # K_CUT_TERM allows for entries with ' + \
                'three or four integers', ofolder, verbosity)
        abort(f, filename)
      try:
        i0 = val[0]
        i1 = int(val[1])
        i2 = int(val[2])
        if len(val) == 3:
          i3 = i2
          i2 = i1
        else:
          i3 = int(val[3])
      except:
        verbose(' # K_CUT_TERM allows for entries with ' + \
                'three or four integers', ofolder, verbosity)
        abort(f, filename)
      # Sanity
      if i1 > i2:
        ia = i2
        i2 = i1
        i1 = ia
      # Label found?
      if i0 not in atom_lab_act:
          verbose(' # Label ' + i0 + ' in K_CUT_TERM not found ' + \
                  'in the list of active atom labels', \
                  ofolder, verbosity)
          abort(f, filename)
      else:
          i0 = atom_lab_act.index(i0) + 1
      ranges.append([i0,i1,i2,i3])
    # Output
    f.write('{0}\n'.format(len(ranges)))
    for rang in ranges:
      f.write('{0} {1} {2} {3}\n'.format(*rang))
    check = 1
  if check == 0:
    f.write('0\n')

  # K_RAD
  check = 0
  if 'K_RAD' in Dictionary:
    val = Dictionary['K_RAD'][0]
    try:
      f.write('{0:7d}\n'.format(int(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:7d}\n'.format(-1))

  # MEMOJ
  check = 0
  if 'MEMOJ' in Dictionary:
    val = Dictionary['MEMOJ'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('Y\n')

  # PIRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'PIRAM' in Dictionary:
      val = Dictionary['PIRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('N\n')

  # VOI_TYPE
  check = 0
  voitypes = {'H82':'0','A67':'1','HAW78':'2','GAUSS':'3'}
  if 'VOI_TYPE' in Dictionary:
    val = Dictionary['VOI_TYPE'][0]
    for ty in list(voitypes.keys()):
        if val.upper() in ty:
            f.write(voitypes[ty]+'\n')
            check = 1
            break
  if check == 0:
    f.write('0\n')

  # VOI_IRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'VOI_IRAM' in Dictionary:
      val = Dictionary['VOI_IRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('N\n')

  # LTE_VOI_IRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'LTE_VOI_IRAM' in Dictionary:
      val = Dictionary['LTE_VOI_IRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('N\n')

  # VOI_PRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'VOI_PRAM' in Dictionary:
      val = Dictionary['VOI_PRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # LTE_VOI_PRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'LTE_VOI_PRAM' in Dictionary:
      val = Dictionary['LTE_VOI_PRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('N\n')

  # RAM_LIM
  check = 0
  if 'RAM_LIM' in Dictionary:
    val = Dictionary['RAM_LIM'][0]
    try:
      if int(val) > 0:
          f.write('{0:15d}\n'.format(int(val)))
          check = 1
      else:
          f.write('-1\n')
          check = 1
    except:
      if val.upper() in 'NONE':
          f.write('-1\n')
          check = 1
      else:
          pass
  if check == 0:
    f.write('{0:15d}\n'.format(-1))

  # RAM_REPORT
  if rmode == 0:
    check = 0
    if 'RAM_REPORT' in Dictionary:
      val = Dictionary['RAM_REPORT'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RAMAN
  check = 0
  if 'RAMAN' in Dictionary:
    val = Dictionary['RAMAN'][0]
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
    elif val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
  if check == 0:
    f.write('Y\n')

  # NO_COH_L_TERM
  check = 0
  if 'NO_COH_L_TERM' in Dictionary:
    val = Dictionary['NO_COH_L_TERM'][0]
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
    elif val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # RED_RESTRICT_HEIGHT
  if rmode == -1 or rmode == 0 or rmode == 1:
    check = 0
    if 'RED_RESTRICT_HEIGHT' in Dictionary:
      val = Dictionary['RED_RESTRICT_HEIGHT'][0]
      try:
        val = interpret(val)
        rang = float(val)
        f.write('Y\n')
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        verbose(' # RED_RESTRICT_HEIGHT wrong format', \
                ofolder, verbosity)
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RED_RESTRICT_TAUC
  if rmode == -1 or rmode == 0 or rmode == 1:
    check = 0
    if 'RED_RESTRICT_TAUC' in Dictionary:
      val = Dictionary['RED_RESTRICT_TAUC'][0]
      try:
        val = interpret(val)
        rang = 10e0**float(val)
        f.write('Y\n')
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        verbose(' # RED_RESTRICT_TAUC wrong format', \
                ofolder, verbosity)
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RED_COHW
  check = 0
  if 'RED_COHW' in Dictionary:
    val = Dictionary['RED_COHW'][0]
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('-1\n')
      dcoh = '-1'
      check = 1
    else:
      try:
        val = interpret(val)
        dw = float(val)
        f.write(val+'\n')
        dcoh = val
        check = 1
      except:
        pass
  if check == 0:
    f.write('-1\n')
    dcoh = '-1'

  # REDI_COHW
  check = 0
  if 'REDI_COHW' in Dictionary:
    val = Dictionary['REDI_COHW'][0]
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('-1\n')
      check = 1
    else:
      try:
        val = interpret(val)
        dw = float(val)
        f.write(val+'\n')
        check = 1
      except:
        pass
  if check == 0:
    f.write(dcoh+'\n')

  # RED_INT_MODE
  check = 0
  if 'RED_INT_MOD' in Dictionary:
    val = Dictionary['RED_INT_MODE'][0]
    if 'LIN' in val:
      f.write('0\n')
      check = 1
    elif 'SPL' in val:
      f.write('1\n')
      check = 1
  if check == 0:
    f.write('1\n')

  # RED_MOD
  check = 0
  if 'RED_MOD' in Dictionary:
    val = Dictionary['RED_MOD'][0]
    if val == 'AA':
      f.write('A\n')
      check = 1
    if val == 'AD':
      f.write('D\n')
      check = 1
  if check == 0:
    f.write('A\n')

  # RED_AAINT
  check = 0
  if 'RED_AAINT' in Dictionary:
    val = Dictionary['RED_AAINT'][0]
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
    elif val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # RED_IRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'RED_IRAM' in Dictionary:
      val = Dictionary['RED_IRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('N\n')

  # RED_PRAM
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'RED_PRAM' in Dictionary:
      val = Dictionary['RED_PRAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # RED_NODE
  check = 0
  if 'RED_NODE' in Dictionary:
    val = Dictionary['RED_NODE'][0]
    try:
      f.write('{0:7d}\n'.format(int(val)))
      check = 1
      redn = int(val)
    except:
      pass
  if check == 0:
    f.write('{0:7d}\n'.format(4))
    redn = 4

  # REDI_NODE
  check = 0
  if 'REDI_NODE' in Dictionary:
    val = Dictionary['REDI_NODE'][0]
    try:
      f.write('{0:7d}\n'.format(int(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:7d}\n'.format(redn))

  # RED_RANG
  check = 0
  if 'RED_RANG' in Dictionary:
    val = Dictionary['RED_RANG'][0]
    try:
      val = interpret(val)
      rang = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(3.5))
    rang = 3.5

  # RED_RESO
  check = 0
  if 'RED_RESO' in Dictionary:
    val = Dictionary['RED_RESO'][0]
    try:
      val = interpret(val)
      reso = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(3.5))
    reso = 3.5

  # RED_NEGL
  check = 0
  if 'RED_NEGL' in Dictionary:
    val = Dictionary['RED_NEGL'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      negl = float(val)
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(10000))
    negl = 10000.

  # RED_VLAR
  check = 0
  if 'RED_VLAR' in Dictionary:
    val = Dictionary['RED_VLAR'][0]
    try:
      val = interpret(val)
      vlar = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(7))
    vlar = 7.

  # RED_FSTP
  check = 0
  if 'RED_FSTP' in Dictionary:
    val = Dictionary['RED_FSTP'][0]
    try:
      val = interpret(val)
      fstp = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(0.5))
    fstp = 0.5

  # RED_MSTP
  check = 0
  if 'RED_MSTP' in Dictionary:
    val = Dictionary['RED_MSTP'][0]
    try:
      val = interpret(val)
      mstp = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(4))
    mstp = 4

  # RED_CORE
  check = 0
  if 'RED_CORE' in Dictionary:
    val = Dictionary['RED_CORE'][0]
    try:
      val = interpret(val)
      core = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(0))
    core = 0

  # RED_RANG_CORE
  check = 0
  if 'RED_RANG_CORE' in Dictionary:
    val = Dictionary['RED_RANG_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(rang))

  # RED_VLAR_CORE
  check = 0
  if 'RED_VLAR_CORE' in Dictionary:
    val = Dictionary['RED_VLAR_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(vlar))

  # RED_FSTP_CORE
  check = 0
  if 'RED_FSTP_CORE' in Dictionary:
    val = Dictionary['RED_FSTP_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(fstp))

  # RED_MSTP_CORE
  check = 0
  if 'RED_MSTP_CORE' in Dictionary:
    val = Dictionary['RED_MSTP_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(mstp))

  # REDI_RANG
  check = 0
  if 'REDI_RANG' in Dictionary:
    val = Dictionary['REDI_RANG'][0]
    try:
      val = interpret(val)
      rangi = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(rang))
    rangi = rang

  # REDI_RESO
  check = 0
  if 'REDI_RESO' in Dictionary:
    val = Dictionary['REDI_RESO'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(reso))

  # REDI_NEGL
  check = 0
  if 'REDI_NEGL' in Dictionary:
    val = Dictionary['REDI_NEGL'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(negl))

  # REDI_VLAR
  check = 0
  if 'REDI_VLAR' in Dictionary:
    val = Dictionary['REDI_VLAR'][0]
    try:
      val = interpret(val)
      vlari = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(vlar))
    vlari = vlar

  # REDI_FSTP
  check = 0
  if 'REDI_FSTP' in Dictionary:
    val = Dictionary['REDI_FSTP'][0]
    try:
      val = interpret(val)
      fstpi = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(fstp))
    fstpi = fstp

  # REDI_MSTP
  check = 0
  if 'REDI_MSTP' in Dictionary:
    val = Dictionary['REDI_MSTP'][0]
    try:
      val = interpret(val)
      mstpi = float(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(mstp))
    mstpi = mstp

  # REDI_CORE
  check = 0
  if 'REDI_CORE' in Dictionary:
    val = Dictionary['REDI_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(core))

  # REDI_RANG_CORE
  check = 0
  if 'REDI_RANG_CORE' in Dictionary:
    val = Dictionary['REDI_RANG_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(rangi))

  # REDI_VLAR_CORE
  check = 0
  if 'REDI_VLAR_CORE' in Dictionary:
    val = Dictionary['REDI_VLAR_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(vlari))

  # REDI_FSTP_CORE
  check = 0
  if 'REDI_FSTP_CORE' in Dictionary:
    val = Dictionary['REDI_FSTP_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(fstpi))

  # REDI_MSTP_CORE
  check = 0
  if 'REDI_MSTP_CORE' in Dictionary:
    val = Dictionary['REDI_MSTP_CORE'][0]
    try:
      val = interpret(val)
      f.write('{0:22.16e}\n'.format(float(val)))
      check = 1
    except:
      pass
  if check == 0:
    f.write('{0:22.16e}\n'.format(mstpi))

  # DOP_WIDTH
  check = 0
  if 'DOP_WIDTH' in Dictionary:
    val = Dictionary['DOP_WIDTH'][0]
    try:
      val = interpret(val)
      num = float(val)
      f.write('NUM\n')
      f.write(val+'\n')
      check = 1
    except ValueError:
      if val == 'MAX' or val == 'MIN':
        if rmode == 0:
          f.write(val+'\n')
          check = 1
        else:
          verbose(' # DOP_WIDTH=MIN/MAX only valid in 1D synthesis', \
                  ofolder, verbosity)
    except:
      verbose(sys.exc_info()[0] + '\n' + sys.exc_info()[1], \
              ofolder, verbosity)
      abort(f, filename)
  if check == 0:
    f.write('NUM\n')
    f.write('2.5d3\n')

  # FORCE_MICRO
  if rmode == 0 or rmode == 1 or rmode == 2:
    check = 0
    if 'FORCE_MICRO' in Dictionary:
      val = Dictionary['FORCE_MICRO'][0]
      try:
        val = interpret(val)
        num = float(val)
        if num > 0:
          f.write(val+'\n')
          check = 1
      except:
        pass
    if check == 0:
      f.write('-1\n')
  else:
    f.write('-1\n')

  # MIN_T
  if rmode == -1 or rmode == 1 or rmode == 2:
    check = 0
    if 'MIN_T' in Dictionary:
      val = Dictionary['MIN_T'][0]
      try:
        val = interpret(val)
        num = float(val)
        if num > 0:
          f.write(val+'\n')
          check = 1
      except:
        pass
    if check == 0:
      f.write('-1\n')
  else:
    f.write('-1\n')

  # MAX_T
  if rmode == -1 or rmode == 1 or rmode == 2:
    check = 0
    if 'MAX_T' in Dictionary:
      val = Dictionary['MAX_T'][0]
      try:
        val = interpret(val)
        num = float(val)
        if num > 0:
          f.write(val+'\n')
          check = 1
      except:
        pass
    if check == 0:
      f.write('-1\n')
  else:
    f.write('-1\n')

  # MAX_V
  if rmode == -1 or rmode == 1 or rmode == 2:
    check = 0
    if 'MAX_V' in Dictionary:
      val = Dictionary['MAX_V'][0]
      try:
        val = interpret(val)
        num = float(val)
        if num > 0:
          f.write(val+'\n')
          check = 1
      except:
        pass
    if check == 0:
      f.write('-1\n')
  else:
    f.write('-1\n')

  # RT_GROUP_N
  if rmode == 1 or rmode == 2 or rmode == -1:
    check = 0
    if 'RT_GROUP_N' in Dictionary:
      val = Dictionary['RT_GROUP_N'][0]
      try:
        num = int(val)
        if num > 0:
          f.write(val+'\n')
          check = 1
      except:
        pass
    if check == 0:
      f.write('1\n')
  else:
    f.write('1\n')

  # UNMAGNETIZED
  if rmode == 1 or rmode == 2:
    check = 0
    if 'UNMAGNETIZED' in Dictionary:
      val = Dictionary['UNMAGNETIZED'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # STATIC
  if rmode == 0 or rmode == 1 or rmode == 2:
    check = 0
    if 'STATIC' in Dictionary:
      val = Dictionary['STATIC'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        static = 'Y\n'
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
        static = 'N\n'
    if check == 0:
      f.write('N\n')
      static = 'N\n'
  else:
    f.write('N\n')
    static = 'N\n'

  # STATIC_INT
  if rmode == 0 or rmode == -1 or rmode == 0:
    check = 0
    if 'STATIC_INT' in Dictionary:
      val = Dictionary['STATIC_INT'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write(static)
  else:
    f.write('N\n')

  # SKIP_DISK
  if rmode == 2:
    check = 0
    if 'SKIP_DISK' in Dictionary:
      val = Dictionary['SKIP_DISK'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # INIT_J_BB
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'INIT_J_BB' in Dictionary:
      val = Dictionary['INIT_J_BB'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('Y\n')

  # ITER_MIN
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MIN' in Dictionary:
      val = Dictionary['ITER_MIN'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
    if check == 0:
      f.write('{0:7d}\n'.format(1))
  else:
    f.write('0\n')

  # ITERI_MIN
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITERI_MIN' in Dictionary:
      val = Dictionary['ITERI_MIN'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
    if check == 0:
      f.write('{0:7d}\n'.format(1))
  else:
    f.write('0\n')

  # ITER_MAX
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MAX' in Dictionary:
      val = Dictionary['ITER_MAX'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
        itmax = int(val)
    if check == 0:
      f.write('{0:7d}\n'.format(500))
      itmax = 500
  else:
    f.write('0\n')

  # ITERI_MAX
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITERI_MAX' in Dictionary:
      val = Dictionary['ITERI_MAX'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
    if check == 0:
      f.write('{0:7d}\n'.format(itmax))
  else:
    f.write('0\n')

  # ITER_2ORD
  if rmode >= -1 and rmode <= 2:
    check = 0
    if 'ITER_2ORD' in Dictionary:
      val = Dictionary['ITER_2ORD'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # ITER_NB
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_NB' in Dictionary:
      val = Dictionary['ITER_NB'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # ITER_MRC_J
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MRC_J' in Dictionary:
      val = Dictionary['ITER_MRC_J'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-4))
  else:
    f.write('0\n')

  # ITER_MRC_I
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MRC_I' in Dictionary:
      val = Dictionary['ITER_MRC_I'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        mrci = float(val)
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-5))
      mrci = 1e-5
  else:
    f.write('0\n')

  # ITERI_MRC_I
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITERI_MRC_I' in Dictionary:
      val = Dictionary['ITERI_MRC_I'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(mrci))
  else:
    f.write('0\n')

  # ITER_MRC_P
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MRC_P' in Dictionary:
      val = Dictionary['ITER_MRC_P'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-3))
  else:
    f.write('0\n')

  # ITER_J
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_J' in Dictionary:
      val = Dictionary['ITER_J'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
        itmax = int(val)
    if check == 0:
      f.write('{0:7d}\n'.format(5))
  else:
    f.write('0\n')

  # ITERI_PRD
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITERI_PRD' in Dictionary:
      val = Dictionary['ITERI_PRD'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
        itmax = int(val)
    if check == 0:
      f.write('{0:7d}\n'.format(4))
  else:
    f.write('0\n')

  # ITERI_MRC_R
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITERI_MRC_R' in Dictionary:
      val = Dictionary['ITERI_MRC_R'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        mrci = float(val)
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-3))
  else:
    f.write('0\n')

  # ITER_PRD
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_PRD' in Dictionary:
      val = Dictionary['ITER_PRD'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
        itmax = int(val)
    if check == 0:
      f.write('{0:7d}\n'.format(1))
  else:
    f.write('0\n')

  # ITER_MRC_R
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MRC_R' in Dictionary:
      val = Dictionary['ITER_MRC_R'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        mrci = float(val)
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-2))
  else:
    f.write('0\n')

  # ITER_MRC_P_R
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ITER_MRC_P_R' in Dictionary:
      val = Dictionary['ITER_MRC_P_R'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        mrci = float(val)
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-1))
  else:
    f.write('0\n')

  # NG_ACC
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'NG_ACC' in Dictionary:
      val = Dictionary['NG_ACC'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        NGACC = 'Y'
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        NGACC = 'N'
        check = 1
    if check == 0:
      f.write('N\n')
      NGACC = 'N'
  else:
    f.write('N\n')

  # NG_ORD
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'NG_ORD' in Dictionary:
      val = Dictionary['NG_ORD'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        NGORD = int(val)
        check = 1
      except:
        NGORD = 3
        pass
    if check == 0:
      NGORD = 3
      f.write('{0:7d}\n'.format(3))
  else:
    f.write('0\n')

  # NG_DELAY
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'NG_DELAY' in Dictionary:
      val = Dictionary['NG_DELAY'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        NGDELAY = int(val)
        check = 1
      except:
        NGDELAY = 20
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(20))
      NGDELAY = 20
  else:
    f.write('2000\n')

  # NGI_ACC
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'NGI_ACC' in Dictionary:
      val = Dictionary['NGI_ACC'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write(NGACC+'\n')
  else:
    f.write('N\n')

  # NGI_ORD
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'NGI_ORD' in Dictionary:
      val = Dictionary['NGI_ORD'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(NGORD))
  else:
    f.write('0\n')

  # NGI_DELAY
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'NGI_DELAY' in Dictionary:
      val = Dictionary['NGI_DELAY'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(NGDELAY))
  else:
    f.write('2000\n')

  # PRD_DELAY
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'PRD_DELAY' in Dictionary:
      val = Dictionary['PRD_DELAY'][0]
      try:
        val = int(val)
        if val > 0:
          f.write('{0:7d}\n'.format(val))
          check = 1
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(1))
  else:
    f.write('0\n')

  # ALI_PHOTO
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ALI_PHOTO' in Dictionary:
      val = Dictionary['ALI_PHOTO'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('Y\n')

  # ALI_DELAY
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ALI_DELAY' in Dictionary:
      val = Dictionary['ALI_DELAY'][0]
      try:
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:7d}\n'.format(0))
  else:
    f.write('0\n')

  # ALI_FORCE
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ALI_FORCE' in Dictionary:
      val = Dictionary['ALI_FORCE'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # ALI_ALLOW_OFF
  if rmode >= -1 and rmode <= 1:
    check = 0
    if 'ALI_ALLOW_OFF' in Dictionary:
      val = Dictionary['ALI_ALLOW_OFF'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')
  else:
    f.write('Y\n')

  # APPEND_MRC
  if rmode == 0:
    check = 0
    if 'APPEND_MRC' in Dictionary:
      val = Dictionary['APPEND_MRC'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # APPENDI_MRC
  if rmode == 0:
    check = 0
    if 'APPENDI_MRC' in Dictionary:
      val = Dictionary['APPENDI_MRC'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # ALLOW_NPHYS_STK
  check = 0
  if 'ALLOW_NPHYS_STK' in Dictionary:
    val = Dictionary['ALLOW_NPHYS_STK'][0]
    try:
      val = int(val)
      check = 1
    except:
      pass
    if(check == 1):
      if(val > 0):
        f.write('{0:7d}\n'.format(val))
      else:
        check = 0
  if check == 0:
    f.write('{0:7d}\n'.format(-1))

  # ALLOW_NPHYS_RHO
  check = 0
  if 'ALLOW_NPHYS_RHO' in Dictionary:
    val = Dictionary['ALLOW_NPHYS_RHO'][0]
    try:
      val = int(val)
      check = 1
    except:
      pass
    if(check == 1):
      if(val > 0):
        f.write('{0:7d}\n'.format(val))
      else:
        check = 0
  if check == 0:
    f.write('{0:7d}\n'.format(-1))

  # ALLOW_NPHYS_POP
  check = 0
  if 'ALLOW_NPHYS_POP' in Dictionary:
    val = Dictionary['ALLOW_NPHYS_POP'][0]
    try:
      val = int(val)
      check = 1
    except:
      pass
    if(check == 1):
      if(val > 0):
        f.write('{0:7d}\n'.format(val))
      else:
        check = 0
  if check == 0:
    f.write('{0:7d}\n'.format(-1))

  # SOLUTION_BOX
  check = 0
  if rmode == 1 or rmode == -1:
    if 'SOLUTION_BOX' in Dictionary:
      val = Dictionary['SOLUTION_BOX']
      if isinstance(val, list):
        if len(val) == 4:
          out = ''
          for v in val:
            try:
              f.write(' {0}\n'.format(int(v)))
            except ValueError:
              verbose(' # SOLUTION_BOX elements must be integers', \
                      '', verbosity)
              abort(f, filename)
            except:
              verbose(' # Problem reading SOLUTION_BOX', \
                      '', verbosity)
              abort(f, filename)
        else:
          verbose(' # SOLUTION_BOX must have 4 elements', \
                  '', verbosity)
          abort(f, filename)
        check = 1
      else:
        verbose(' # SOLUTION_BOX must be a list', '', verbosity)
        abort(f, filename)
  if check == 0:
      f.write('-1\n-1\n-1\n-1\n')

  # EXCLUDE_PIXEL
  check = 0
  if rmode == -1 or rmode == 1:
    if 'EXCLUDE_PIXEL' in Dictionary:
      vals = Dictionary['EXCLUDE_PIXEL']
      NL = len(vals)
      pixels = []
      for val in vals:
        if len(val) != 2:
          verbose(' # EXCLUDE_PIXEL allows for entries with ' + \
                  'two integers', ofolder, verbosity)
          abort(f, filename)
        try:
          i0 = int(val[0])
          i1 = int(val[1])
        except:
          verbose(' # EXCLUDE_PIXEL allows for entries with ' + \
                  'two integers', ofolder, verbosity)
          abort(f, filename)
        pixels.append([i0,i1])
      # Process
      pixels = process_pixels(pixels)
      # Output
      f.write('{0}\n'.format(len(pixels)))
      for pix in pixels:
        f.write('{0} {1}\n'.format(*pix))
      check = 1
  if check == 0:
    f.write('0\n')

  # STORE_STEP
  if rmode == 0 or rmode == 1:
    check = 0
    valid = 0
    if 'STORE_STEP' in Dictionary:
      check = 1
      val = Dictionary['STORE_STEP'][0]
      if val.isdigit():
        val = int(val)
        if val > 0:
          f.write('{0:7d}\n'.format(val))
          valid = 1
    if valid == 0 and check == 1:
      f.write('{0:7d}\n'.format(-1))
    if check == 0:
      f.write('{0:7d}\n'.format(-1))
  else:
    f.write('{0:7d}\n'.format(-1))

  # STOREI_STEP
  if rmode == 0 or rmode == 1:
    check = 0
    valid = 0
    if 'STOREI_STEP' in Dictionary:
      check = 1
      val = Dictionary['STOREI_STEP'][0]
      if val.isdigit():
        val = int(val)
        if val > 0:
          f.write('{0:7d}\n'.format(val))
          valid = 1
    if valid == 0 and check == 1:
      f.write('{0:7d}\n'.format(-1))
    if check == 0:
      f.write('{0:7d}\n'.format(-1))
  else:
    f.write('{0:7d}\n'.format(-1))

  # CONTRIBUTION
  keep_ctr = False
  if rmode == -1 or rmode == 0 or rmode == 1:
    check = 0
    if 'CONTRIBUTION' in Dictionary:
      val = Dictionary['CONTRIBUTION'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        keep_ctr = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # TAU1
  keep_tau = False
  if rmode >= -1:
    check = 0
    if 'TAU1' in Dictionary:
      val = Dictionary['TAU1'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        keep_tau = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_BACK
  keep_back = False
  if rmode >= 0:
    check = 0
    if 'KEEP_BACK' in Dictionary:
      val = Dictionary['KEEP_BACK'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        keep_back = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_DAMP
  keep_damp = False
  if rmode >= 0:
    check = 0
    if 'KEEP_DAMP' in Dictionary:
      val = Dictionary['KEEP_DAMP'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        keep_damp = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_QEL
  keep_qel = False
  if rmode >= 0:
    check = 0
    if 'KEEP_QEL' in Dictionary:
      val = Dictionary['KEEP_QEL'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        keep_qel = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_APARAM
  if rmode == 0:
    check = 0
    if 'KEEP_APARAM' in Dictionary:
      val = Dictionary['KEEP_APARAM'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_COLS
  keep_cols = False
  if rmode >= 0:
    check = 0
    if 'KEEP_COLS' in Dictionary:
      val = Dictionary['KEEP_COLS'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        keep_cols = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_ATMO
  if rmode >= 0:
    check = 0
    if 'KEEP_ATMO' in Dictionary:
      val = Dictionary['KEEP_ATMO'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_POP
  if rmode >= 0:
    check = 0
    if 'KEEP_POP' in Dictionary:
      val = Dictionary['KEEP_POP'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      if rmode == 0:
        f.write('Y\n')
      else:
        f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_DEP
  if rmode >= 0:
    check = 0
    if 'KEEP_DEP' in Dictionary:
      val = Dictionary['KEEP_DEP'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      if rmode == 0:
        f.write('Y\n')
      else:
        f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_RHOKQ
  if rmode >= 0:
    check = 0
    if 'KEEP_RHOKQ' in Dictionary:
      val = Dictionary['KEEP_RHOKQ'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      if rmode == 0:
        f.write('Y\n')
      else:
        f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_JKQ
  if rmode >= 0:
    check = 0
    if 'KEEP_JKQ' in Dictionary:
      val = Dictionary['KEEP_JKQ'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      if rmode == 0:
        f.write('Y\n')
      else:
        f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_STOKES_QUAD
  if rmode >= 0:
    check = 0
    if 'KEEP_STOKES_QUAD' in Dictionary:
      val = Dictionary['KEEP_STOKES_QUAD'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      if rmode == 0:
        f.write('Y\n')
      else:
        f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_JKQNU
  if rmode >= 0:
    check = 0
    if 'KEEP_JKQNU' in Dictionary:
      val = Dictionary['KEEP_JKQNU'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_MRC
  if rmode == 0 or rmode == 1:
    check = 0
    if 'KEEP_MRC' in Dictionary:
      val = Dictionary['KEEP_MRC'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_COL_LOG
  if rmode == 0:
    check = 0
    if 'KEEP_COL_LOG' in Dictionary:
      val = Dictionary['KEEP_COL_LOG'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # KEEP_MPI_LOG
  check = 0
  if 'KEEP_MPI_LOG' in Dictionary:
    val = Dictionary['KEEP_MPI_LOG'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # KEEP_MPI_DETAIL_LOG
  check = 0
  if 'KEEP_MPI_DETAIL_LOG' in Dictionary:
    val = Dictionary['KEEP_MPI_DETAIL_LOG'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # LIM_STK
  if rmode == 1 or rmode == 2:
    if 'LIM_STK' in Dictionary:
      vals = Dictionary['LIM_STK']
      NL = len(vals)
      doublets = []
      for val in vals:
        try:
          val[0] = interpret(val[0])
          val[1] = interpret(val[1])
          l0 = 1e2/float(val[0])
          l1 = 1e2/float(val[1])
          if l1 >= l0:
            doublets.append([l0,l1])
          else:
            doublets.append([l1,l0])
        except:
          verbose(' # LIM_STK wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      NL, doublets = Worder(NL,doublets)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_CTR
  if rmode == 1 and keep_ctr:
    if 'LIM_CTR' in Dictionary:
      vals = Dictionary['LIM_CTR']
      NL = len(vals)
      doublets = []
      for val in vals:
        try:
          val[0] = interpret(val[0])
          val[1] = interpret(val[1])
          l0 = 1e2/float(val[0])
          l1 = 1e2/float(val[1])
          if l1 >= l0:
            doublets.append([l0,l1])
          else:
            doublets.append([l1,l0])
        except:
          verbose(' # LIM_CTR wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      NL, doublets = Worder(NL,doublets)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_TAU
  if rmode == 1 and keep_tau:
    if 'LIM_TAU' in Dictionary:
      vals = Dictionary['LIM_TAU']
      NL = len(vals)
      doublets = []
      for val in vals:
        try:
          val[0] = interpret(val[0])
          val[1] = interpret(val[1])
          l0 = 1e2/float(val[0])
          l1 = 1e2/float(val[1])
          if l1 >= l0:
            doublets.append([l0,l1])
          else:
            doublets.append([l1,l0])
        except:
          verbose(' # LIM_TAU wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      NL, doublets = Worder(NL,doublets)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_COLS_TT
  if rmode == 1 and keep_cols:
    if 'LIM_COLS_TT' in Dictionary:
      vals = Dictionary['LIM_COLS_TT']
      NL = len(vals)
      triplets = []
      for val in vals:
        ia = val[0]
        # Label found?
        if ia not in atom_lab_act:
            verbose(' # Label ' + ia + ' in LIM_COLS_TT not ' + \
                    'found in the list of active atom labels', \
                    ofolder, verbosity)
            abort(f, filename)
        else:
            ia = atom_lab_act.index(ia) + 1
        try:
          i0 = int(val[1])
          i1 = int(val[2])
          if i0 > i1:
            triplets.append([ia,i0,i1])
          else:
            triplets.append([ia,i1,i0])
        except:
          verbose(' # LIM_COLS_TT wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      f.write('{0}\n'.format(NL))
      for tri in triplets:
        f.write('{0}\n'.format(tri[0]))
        f.write('{0}\n'.format(tri[1]))
        f.write('{0}\n'.format(tri[2]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_COLS_LL
  if rmode == 1 and keep_cols:
    if 'LIM_COLS_LL' in Dictionary:
      vals = Dictionary['LIM_COLS_LL']
      NL = len(vals)
      triplets = []
      for val in vals:
        ia = val[0]
        # Label found?
        if ia not in atom_lab_act:
            verbose(' # Label ' + ia + ' in LIM_COLS_LL not ' + \
                    'found in the list of active atom labels', \
                    ofolder, verbosity)
            abort(f, filename)
        else:
            ia = atom_lab_act.index(ia) + 1
        try:
          i0 = int(val[1])
          i1 = int(val[2])
          if i0 > i1:
            triplets.append([ia,i0,i1])
          else:
            triplets.append([ia,i1,i0])
        except:
          verbose(' # LIM_COLS_LL wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      f.write('{0}\n'.format(NL))
      for tri in triplets:
        f.write('{0}\n'.format(tri[0]))
        f.write('{0}\n'.format(tri[1]))
        f.write('{0}\n'.format(tri[2]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_DAMP
  if rmode == 1 and keep_damp:
    if 'LIM_DAMP' in Dictionary:
      vals = Dictionary['LIM_DAMP']
      NL = len(vals)
      doublets = []
      for val in vals:
        ia = val[0]
        # Label found?
        if ia not in atom_lab_act:
            verbose(' # Label ' + ia + ' in LIM_DAMP not ' + \
                    'found in the list of active atom labels', \
                    ofolder, verbosity)
            abort(f, filename)
        else:
            ia = atom_lab_act.index(ia) + 1
        try:
          it = int(val[1])
          doublets.append([ia,it])
        except:
          verbose(' # LIM_DAMP wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_QEL
  if rmode == 1 and keep_qel:
    if 'LIM_QEL' in Dictionary:
      vals = Dictionary['LIM_QEL']
      NL = len(vals)
      doublets = []
      for val in vals:
        ia = val[0]
        # Label found?
        if ia not in atom_lab_act:
            verbose(' # Label ' + ia + ' in LIM_QEL not ' + \
                    'found in the list of active atom labels', \
                    ofolder, verbosity)
            abort(f, filename)
        else:
            ia = atom_lab_act.index(ia) + 1
        try:
          it = int(val[1])
          doublets.append([ia,it])
        except:
          verbose(' # LIM_QEL wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_BACK
  if rmode == 1 and keep_back:
    if 'LIM_BACK' in Dictionary:
      vals = Dictionary['LIM_BACK']
      NL = len(vals)
      doublets = []
      for val in vals:
        try:
          val[0] = interpret(val[0])
          val[1] = interpret(val[1])
          l0 = 1e2/float(val[0])
          l1 = 1e2/float(val[1])
          if l1 >= l0:
            doublets.append([l0,l1])
          else:
            doublets.append([l1,l0])
        except:
          verbose(' # LIM_BACK wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # LIM_POP
  if rmode == 1:
    if 'LIM_POP' in Dictionary:
      vals = Dictionary['LIM_POP']
      NL = len(vals)
      doublets = []
      for val in vals:
        ia = val[0]
        # Label found?
        if ia not in atom_lab_act:
            verbose(' # Label ' + ia + ' in LIM_POP not ' + \
                    'found in the list of active atom labels', \
                    ofolder, verbosity)
            abort(f, filename)
        else:
            ia = atom_lab_act.index(ia) + 1
        try:
          it = int(val[1])
          doublets.append([ia,it])
        except:
          verbose(' # LIM_POP wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      f.write('{0}\n'.format(NL))
      for dou in doublets:
        f.write('{0}\n'.format(dou[0]))
        f.write('{0}\n'.format(dou[1]))
    else:
      f.write('0\n')
  else:
    f.write('0\n')

  # REDO_NE
  check = 0
  if 'REDO_NE' in Dictionary:
    val = Dictionary['REDO_NE'][0]
    if 'Y' in val or 'S' in val:
      f.write('10\n')
      check = 10
    if val in 'INIT':
      f.write('10\n')
      check = 1
    if val in 'FIN':
      f.write('1\n')
      check = 1
    if val in 'BOTH':
      f.write('11\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('0\n')
      check = 1
  if check == 0:
    f.write('0\n')

  # UPDATE_ATMOS
  check = 0
  if rmode == 0:
    if 'UPDATE_ATMOS' in Dictionary:
      val = Dictionary['UPDATE_ATMOS'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('0\n')
        check = 1
      elif val == 'S' or val == 'ST' or val == 'STA' or \
           val == 'STAN' or val == 'STAND' or \
           val == 'STANDA' or val == 'STANDARD':
        f.write('1\n')
        check = 1
      elif val == 'N' or val == 'NE':
        f.write('2\n')
        check = 1
      elif val == 'PE':
        f.write('3\n')
        check = 1
      elif val == 'RE' or val == 'RHOE':
        f.write('4\n')
        check = 1
      elif val == 'P' or val == 'PG' or val == 'PGA' or val == 'PGAS':
        f.write('5\n')
        check = 1
      elif val == 'R' or val == 'RH' or val == 'RHO':
        f.write('6\n')
        check = 1
      else:
        f.write('-1\n')
        check = 1
    if check == 0:
      f.write('-1\n')
  else:
    f.write('-1\n')

  # PROTECT_H
  check = 0
  if 'PROTECT_H' in Dictionary:
    val = Dictionary['PROTECT_H'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # CHEM_PROTECT_ALL
  check = 0
  if 'CHEM_PROTECT_ALL' in Dictionary:
    val = Dictionary['CHEM_PROTECT_ALL'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  # WRITE_PERFORMANCE
  if rmode == 0:
    check = 0
    if 'WRITE_PERFORMANCE' in Dictionary:
      val = Dictionary['WRITE_PERFORMANCE'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # WRITE_MPI_PERFORMANCE
  if rmode == 0:
    check = 0
    if 'WRITE_MPI_PERFORMANCE' in Dictionary:
      val = Dictionary['WRITE_MPI_PERFORMANCE'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')
  else:
    f.write('N\n')

  # VERBOSE
  check = 0
  if 'VERBOSE' in Dictionary:
    val = Dictionary['VERBOSE'][0]
    if val == 'Y' or val == 'YE' or val == 'YES' or \
       val == 'S' or val =='SI':
      f.write('Y\n')
      check = 1
    if val == 'N' or val == 'NO' or val == 'NON':
      f.write('N\n')
      check = 1
  if check == 0:
    f.write('N\n')

  ###################################################################
  # TIC - inversion only keywords
  ###################################################################

  if rmode == -1:

    # VERBOSE_INV_LV
    check = 0
    if 'VERBOSE_INV_LV' in Dictionary:
      allowed = [0,1,2,3]
      val = Dictionary['VERBOSE_INV_LV'][0]
      try:
        if int(val) in allowed:
          f.write(val+'\n')
          check = 1
      except ValueError:
        verbose('VERBOSE_INV_LV must be an integer', \
                ofolder,verbosity)
      except:
        pass
    if check == 0:
      f.write('0\n')

    # VERBOSE_INV_SHUTUP
    check = 0
    if 'VERBOSE_INV_SHUTUP' in Dictionary:
      allowed = [0,1,2,3]
      val = Dictionary['VERBOSE_INV_SHUTUP'][0]
      try:
        if int(val) in allowed:
          f.write(val+'\n')
          check = 1
      except ValueError:
        verbose('VERBOSE_INV_SHUTUP must be an integer', \
                ofolder,verbosity)
      except:
        pass
    if check == 0:
      f.write('3\n')

    # DATA_FILE
    # Existence was checked above
    val = Dictionary['DATA_FILE'][0]
    f.write(val+'\n')

    # TYPE_INVERSION
    check = 0
    if 'TYPE_INVERSION' in Dictionary:
      val = Dictionary['TYPE_INVERSION'][0]
      if 'T' in val and 'ET' not in val and 'S' not in val:
        f.write('0\n')
        check = 1
      elif 'B' in val or 'M' in val and 'S' not in val:
        f.write('1\n')
        check = 1
      elif 'A' in val and 'S' not in val:
        f.write('2\n')
        check = 1
      elif 'S' in val and 'M' in val:
        f.write('4\n')
        check = 1
      elif 'S' in val:
        f.write('3\n')
        check = 1
    if check == 0:
      f.write('0\n')

    # AUTO_WEIGHT
    check = 0
    if 'AUTO_WEIGHT' in Dictionary:
      val = Dictionary['AUTO_WEIGHT'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
        aweig = True
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
        aweig = False
    if check == 0:
      f.write('N\n')
      aweig = False

    # If weight not automatic
    if not aweig:

      # Initialize bool to check if any weight
      anyweight = False

      # WEIGHT
      check = 0
      if 'WEIGHT' in Dictionary:
        vals = Dictionary['WEIGHT']
        # If only one
        if len(vals) == 1:
          # If size 1
          if len(vals[0]) == 4:
            try:
              val1 = interpret(vals[0][0])
              num = float(val1)
              val2 = interpret(vals[0][1])
              num = float(val2)
              val3 = interpret(vals[0][2])
              num = float(val3)
              val4 = interpret(vals[0][3])
              num = float(val4)
              f.write('1\n'+val1+'\n'+ \
                            val2+'\n'+ \
                            val3+'\n'+ \
                            val4+'\n'+ \
                            '0\n1d100\n')
              check = 1
              anyweight = True
            except:
              verbose(' # WEIGHT must be four or ' + \
                      'six floats', ofolder, verbosity)
              abort(f, filename)
        # Not global weight
        if check == 0:
          NL = len(vals)
          low = []
          up = []
          ff = []
          for val in vals:
            try:
              val[0] = interpret(val[0])
              val[1] = interpret(val[1])
              val[2] = interpret(val[2])
              val[3] = interpret(val[3])
              val[4] = interpret(val[4])
              val[5] = interpret(val[5])
              l = float(val[0])
              u = float(val[1])
              c1 = float(val[2])
              c2 = float(val[3])
              c3 = float(val[4])
              c4 = float(val[5])
              if l > u:
                low.append(u)
                up.append(l)
                ff.append([c1,c2,c3,c4])
              elif u > l:
                low.append(l)
                up.append(u)
                ff.append([c1,c2,c3,c4])
              else:
                verbose(' # WEIGHT ranges must have size ' + \
                        'larger than 0 nm', ofolder, verbosity)
                abort(f, filename)
            except:
              verbose(' # WEIGHT must be four or six floats', \
                      ofolder, verbosity)
              abort(f, filename)
          valid,n,low,up,ff = unique_ranges(NL,low,up,ff,[0.,1e100])
          if valid:
            if n > 0:
              f.write('{0}\n'.format(n))
              msg = '{0} {1} {2} {3} {4} {5}\n'
              for l,u,c in zip(low,up,ff):
                  f.write(msg.format(c[0],c[1],c[2],c[3],l,u))
              check = 1
              anyweight = True
            else:
              f.write('0\n')
              check = 1
          else:
            verbose(' # WEIGHT ranges must not intersect', \
                    ofolder, verbosity)
            abort(f, filename)
      if check == 0:
        f.write('0\n')

      # WEIGHT_FILE
      if 'WEIGHT_FILE' in Dictionary:
        val = Dictionary['WEIGHT_FILE'][0]
        try:
          f2=open(val)
          f2.close()
          f.write('1\n'+val+'\n')
          anyweight = True
        except:
          f.write('0\n')
      else:
        f.write('0\n')

      # If no weight defined
      if not anyweight:
        verbose(' # WEIGHT or WEIGHT_FILE are necessary ' + \
                'keywords if AUTO_WEIGHT = No', ofolder, verbosity)
        abort(f, filename)

      # WEIGHT_FACTOR
      conversion = {'I': '0.','Q': '1.','U': '2.','V': '3.'}
      check = 0
      if 'WEIGHT_FACTOR' in Dictionary:
        val = Dictionary['WEIGHT_FACTOR']
        out = []
        for lval in val:
          if lval[0].strip().lower() == 'none': continue
          if len(lval) == 4:
            lval[0] = lval[0].strip().upper()
            if lval[0] not in conversion:
              verbose(' # Identifier in WEIGHT_FACTOR ' + \
                      lval[0] + 'not in [I,Q,U,V]', \
                      ofolder, verbosity)
              abort(f, filename)
            lval[0] = conversion[lval[0]]
            try:
              l0 = float(lval[1])
              l1 = float(lval[2])
              ff = float(lval[3])
            except ValueError:
              verbose(' # Arguments in WEIGHT_FACTOR must be ' + \
                      'one string and three floats', \
                      ofolder, verbosity)
              abort(f, filename)
            except:
              raise
            if ff < 0:
              verbose(' # Factor in WEIGHT_FACTOR must be ' + \
                      'non-negative', \
                      ofolder, verbosity)
              abort(f, filename)
            if l0 > l1:
              out.append('{0} {1} {2} {3}\n'.format(lval[0],l1,l0,ff))
            else:
              out.append('{0} {1} {2} {3}\n'.format(lval[0],l0,l1,ff))
          else:
            verbose(' # Number of elements in WEIGHT_FACTOR ' + \
                    'not equal to 4 in line: ' + ' '.join(lval), \
                    ofolder, verbosity)
            abort(f, filename)
        check = 1
        f.write('{0}\n'.format(len(out)))
        if len(out) > 0:
          for ou in out:
            f.write(ou)
      if check == 0:
        f.write('0\n')

    # SIGMA_FACTOR
    conversion = {'I': '0.','Q': '1.','U': '2.','V': '3.'}
    check = 0
    if 'SIGMA_FACTOR' in Dictionary:
      val = Dictionary['SIGMA_FACTOR']
      out = []
      for lval in val:
        if lval[0].strip().lower() == 'none': continue
        if len(lval) == 4:
          lval[0] = lval[0].strip().upper()
          if lval[0] not in conversion:
            verbose(' # Identifier in SIGMA_FACTOR ' + \
                    lval[0] + 'not in [I,Q,U,V]', \
                    ofolder, verbosity)
            abort(f, filename)
          lval[0] = conversion[lval[0]]
          try:
            l0 = float(lval[1])
            l1 = float(lval[2])
            ff = float(lval[3])
          except ValueError:
            verbose(' # Arguments in SIGMA_FACTOR must be ' + \
                    'one string and three floats', \
                    ofolder, verbosity)
            abort(f, filename)
          except:
            raise
          if ff < 0:
            verbose(' # Factor in SIGMA_FACTOR must be ' + \
                    'non-negative', \
                    ofolder, verbosity)
            abort(f, filename)
          if l0 > l1:
            out.append('{0} {1} {2} {3}\n'.format(lval[0],l1,l0,ff))
          else:
            out.append('{0} {1} {2} {3}\n'.format(lval[0],l0,l1,ff))
        else:
          verbose(' # Number of elements in SIGMA_FACTOR ' + \
                  'not equal to 4 in line: ' + ' '.join(lval), \
                  ofolder, verbosity)
          abort(f, filename)
      check = 1
      f.write('{0}\n'.format(len(out)))
      if len(out) > 0:
        for ou in out:
          f.write(ou)
    if check == 0:
      f.write('0\n')

    # INV_INIT
    check = 0
    if 'INV_INIT' in Dictionary:
      val = Dictionary['INV_INIT'][0]
      if val.upper() in 'INITIALIZE':
        f.write('INIT\n')
        check = 1
      else:
        try:
          f2=open(val)
          f2.close()
        except:
          verbose(' # INV_INIT file not found '+val, \
                  ofolder, verbosity)
          abort(f, filename)
        f.write(val+'\n')
        check = 1
    if check == 0:
      f.write('INIT\n')

    # INV_MASK
    check = 0
    if 'INV_MASK' in Dictionary:
      val = Dictionary['INV_MASK'][0]
      if val.upper() in 'NONE':
        f.write('NONE\n')
        check = 1
      else:
        try:
          f2=open(val)
          f2.close()
        except:
          verbose(' # INV_MASK file not found '+val, \
                  ofolder, verbosity)
          abort(f, filename)
        f.write(val+'\n')
        check = 1
    if check == 0:
      f.write('NONE\n')

    # CENTERED_DERIVATIVE
    check = 0
    if 'CENTERED_DERIVATIVE' in Dictionary:
      val = Dictionary['CENTERED_DERIVATIVE'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

    # ITER_MAX_INV
    check = 0
    if 'ITER_MAX_INV' in Dictionary:
      val = Dictionary['ITER_MAX_INV'][0]
      if val.isdigit():
        f.write('{0:7d}\n'.format(int(val)))
        check = 1
    if check == 0:
      f.write('{0:7d}\n'.format(15))

    # NODES_X_METHOD
    # For each variable
    for var in varis:
      check = 0
      if 'NODES_'+var+'_METHOD' in Dictionary:
        # Special diffuse light
        if var == 'F':
            f.write('0\n')
            continue
        val = ' '.join(Dictionary['NODES_'+var+'_METHOD'])
        if 'V' in val:
          if 'FIX' in val:
            if 'FIR' in val:
              f.write('2\n')
              check = 1
            elif 'LAS' in val:
              f.write('1\n')
              check = 1
            elif 'EXT' in val:
              f.write('3\n')
              check = 1
          else:
            f.write('0\n')
            check = 1
        elif 'C' in val:
          if 'FIX' in val:
            if 'FIR' in val:
              f.write('6\n')
              check = 1
            elif 'LAS' in val:
              f.write('5\n')
              check = 1
            elif 'EXT' in val:
              f.write('7\n')
              check = 1
          else:
            f.write('4\n')
            check = 1
        else:
          msg = ' # NODES_'+var+'_METHOD must be '
          msg += 'value or correction, and additionaly fix ' + \
                 'first, fix last, or fix extremes'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        f.write('0\n')

    # INTERPOLATION
    check = 0
    if 'INTERPOLATION' in Dictionary:
      val = Dictionary['INTERPOLATION'][0]
      if val in 'LINEAR':
        f.write('0\n')
        check = 1
     #elif val in 'QUADRATIC' and val in 'BEZIER':
      elif 'Q' in val and 'B' in val:
        f.write('1\n')
        check = 1
     #elif val in 'CUBIC' and val in 'BEZIER':
      elif 'C' in val and 'B' in val:
        f.write('2\n')
        check = 1
      else:
        msg = ' # INTERPOLATION must be '
        msg += 'linear, quadratic bezier, or cubic bezier'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('2\n')

    # NODES_X_LOCATION/NUM
    # For each variable
    for var in varis:
      check = 0
      if 'NODES_'+var+'_LOCATION' in Dictionary:
        if 'NODES_'+var+'_NUM' in Dictionary:
          msg = ' # NODES_'+var+'_LOCATION and '
          msg += 'NODES_'+var+'_NUM cannot be '
          msg += 'set simultaneously'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
        # Special diffuse light
        if var == 'F':
          msg = ' # You specified NODES_'+var+'_LOCATION, '
          msg += 'setting NODES_'+var+'_NUM = 1'
          verbose(msg, ofolder, verbosity)
          f.write('0\n1\n')
          continue
        vals = Dictionary['NODES_'+var+'_LOCATION']
        nums = order_float(vals)
        if len(nums) > 0:
          f.write('1\n')
          f.write('{0}\n'.format(len(nums)))
          for num in nums:
            f.write(' {0}'.format(num))
            if num > 2:
              msg = ' # Warning: Node location in ' + \
                    var + ' at tau > 2'
              verbose(msg, ofolder, verbosity)
            if num < -9:
              msg = ' # Warning: Node location in ' + \
                    var + ' at tau < -9'
              verbose(msg, ofolder, verbosity)
          f.write('\n')
          check = 1
        else:
          msg = ' # NODES_'+var+'_LOCATION must be '
          msg += 'a list of unique floats'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      elif 'NODES_'+var+'_NUM' in Dictionary:
        val = Dictionary['NODES_'+var+'_NUM'][0]
        try:
          num = int(val)
          if num < 0:
            msg = ' # NODES_'+var+'_NUM must be '
            msg += 'a positive integer'
            verbose(msg, ofolder, verbosity)
            abort(f, filename)
          # Special diffuse light
          if var == 'F':
            if num > 1:
              msg = ' # Diffuse light only admits one node, '
              msg += 'setting NODES_'+var+'_NUM = 1'
              verbose(msg, ofolder, verbosity)
              f.write('0\n1\n')
              continue
          f.write('0\n')
          f.write(val+'\n')
          check = 1
        except:
          msg = ' # NODES_'+var+'_NUM must be '
          msg += 'a positive integer'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        if var not in varis_hidden:
          verbose(' # No node information for variable ' + \
                  var, ofolder, verbosity)
        f.write('0\n')
        f.write('0\n')

    # BTYPE
    check = 0
    if 'BTYPE' in Dictionary:
      val = Dictionary['BTYPE'][0]
      if val in 'VERTICAL':
        f.write('0\n')
        check = 1
        btyp = 0
      elif val in 'LOS' or (val in 'LINE' and val in 'SIGHT'):
        f.write('1\n')
        check = 1
        btyp = 1
      else:
        msg = ' # BTYPE must be vertical or LOS'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0\n')
      btyp = 0

    # VTYPE
    check = 0
    if 'VTYPE' in Dictionary:
      val = Dictionary['VTYPE'][0]
      if val in 'VERTICAL':
        f.write('0\n')
        check = 1
        vtyp = 0
      elif val in 'LOS' or (val in 'LINE' and val in 'SIGHT'):
        f.write('1\n')
        check = 1
        vtyp = 1
      else:
        msg = ' # VTYPE must be vertical or LOS'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0\n')
      vtyp = 0

    # FIX_X
    # For each variable
    for var in varis:
      check = 0
      if 'FIX_'+var in Dictionary:
        val = Dictionary['FIX_'+var][0]
        if val == 'Y' or val == 'YE' or val == 'YES' or \
           val == 'S' or val =='SI':
          f.write('Y\n')
          check = 1
        if val == 'N' or val == 'NO' or val == 'NON':
          f.write('N\n')
          check = 1
      if check == 0:
        f.write('N\n')

    # POSITION_CORRECTION
    check = 0
    if 'POSITION_CORRECTION' in Dictionary:
      val = Dictionary['POSITION_CORRECTION'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')

    # REGUL_X
    # For each variable
    for var in varis:
      check = 0
      oumsg = True
      if 'REGUL_'+var in Dictionary:
        val = Dictionary['REGUL_'+var]
        try:
          if val[0] in 'NONE':
            f.write('0\n')
            check = 1
          elif val[0] in 'MEAN':
            val[1] = interpret(val[1])
            weig = float(val[1])
            f.write('1\n')
            f.write('{0}\n'.format(weig))
            check = 1
          elif val[0] in 'CONSTANT':
            val[1] = interpret(val[1])
            weig = float(val[1])
            f.write('2\n')
            f.write('{0}\n'.format(weig))
            check = 1
          elif 'F' in val[0] and 'D' in val[0]:
            val[1] = interpret(val[1])
            weig = float(val[1])
            f.write('3\n')
            f.write('{0}\n'.format(weig))
            check = 1
          elif 'S' in val[0] and 'D' in val[0]:
            val[1] = interpret(val[1])
            weig = float(val[1])
            f.write('4\n')
            f.write('{0}\n'.format(weig))
            check = 1
          elif 'C' in val[0] and 'L1' in val[0]:
            val[1] = interpret(val[1])
            weig = float(val[1])
            f.write('5\n')
            f.write('{0}\n'.format(weig))
            check = 1
          else:
            msg = ' # REGUL_'+var+' must be '
            msg += 'none, mean, constant, first derivative, '
            msg += 'second derivative, or constantl1'
            verbose(msg, ofolder, verbosity)
            oumsg = False
            abort(f, filename)
        except:
          if oumsg:
            msg = ' # REGUL_'+var+' must be '
            msg += 'none, mean, constant, first derivative, or '
            msg += 'second derivative, followed by the weight if '
            msg += 'not none'
            verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        f.write('0\n')

    # REGUL_LIMITS
    check = 0
    if 'REGUL_LIMITS' in Dictionary:
      val = Dictionary['REGUL_LIMITS'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # REGUL_LIMITS must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0.1\n')

    # REGUL_FACTOR
    check = 0
    if 'REGUL_FACTOR' in Dictionary:
      val = Dictionary['REGUL_FACTOR'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # REGUL_FACTOR must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('1.0\n')

    # THRH_CHI2
    check = 0
    if 'THRH_CHI2' in Dictionary:
      val = Dictionary['THRH_CHI2'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-4))

    # INV_MRC
    check = 0
    if 'INV_MRC' in Dictionary:
      val = Dictionary['INV_MRC'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-2))

    # SVD_TYPE
    check = 0
    if 'SVD_TYPE' in Dictionary:
      val = Dictionary['SVD_TYPE'][0]
      if val in 'TRADITIONAL':
        f.write('0\n')
        check = 1
      elif val in 'SIR':
        f.write('2\n')
        check = 1
      else:
        verbose(' # SVD_TYPE must be traiditional or sir', \
                ofolder, verbosity)
    if check == 0:
      f.write('2\n')

    # THRH_SVD
    check = 0
    if 'THRH_SVD' in Dictionary:
      val = Dictionary['THRH_SVD'][0]
      try:
        val = interpret(val)
        f.write('{0:22.16e}\n'.format(float(val)))
        check = 1
      except:
        pass
    if check == 0:
      f.write('{0:22.16e}\n'.format(1e-4))

    # PG_TYPE
    check = 0
    if 'PG_TYPE' in Dictionary:
      val = Dictionary['PG_TYPE'][0]
     #if val in 'HYDROSTATIC' or val in 'EQUILIBRIUM':
      if 'HY' in val or 'EQ' in val:
        f.write('1\n')
        check = 1
      else:
        f.write('0\n')
        check = 1
    if check == 0:
      f.write('1\n')

    # PG_BOUND
    check = 0
    if 'PG_BOUND' in Dictionary:
      val = Dictionary['PG_BOUND'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # PG_BOUND must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('-1d0\n')

    # DIFFUSE_LIGHT
    check = 0
    if 'DIFFUSE_LIGHT' in Dictionary:
      val = Dictionary['DIFFUSE_LIGHT'][0]
      try:
        val = interpret(val)
        num = float(val)
        if num < 0. or num > 1.:
          verbose(' # DIFFUSE LIGHT must be a float between ' + \
                  '0 and 1', ofolder, verbosity)
          abort(f, filename)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # DIFFUSE LIGHT must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0.\n')

    # ATMO_INPUT
    if 'ATMO_INPUT' in Dictionary:
      val = Dictionary['ATMO_INPUT'][0]
      if val.upper() in 'NONE':
        f.write('NONE\n')
      elif val.upper() == 'HCC':
        f.write('$$C$$\n')
      elif val.upper() == 'HCP':
        f.write('$$P$$\n')
      else:
        try:
          f2=open(val)
          f2.close()
        except:
          verbose(' # ATMO_INPUT file not found '+val, \
                  ofolder, verbosity)
          abort(f, filename)
        f.write(val+'\n')
    else:
      f.write('NONE\n')

    # ATMO_NODES
    check = 0
    if 'ATMO_NODES' in Dictionary:
      val = Dictionary['ATMO_NODES'][0]
      try:
        num = int(val)
        if num < 0:
          verbose('ATMO_NODES must be a positive integer', \
                  ofolder,verbosity)
        f.write(val+'\n')
        check = 1
      except:
        verbose('ATMO_NODES must be a positive integer', \
                ofolder,verbosity)
    if check == 0:
      f.write('0\n')

    # MAX_SVD_STEP
    check = 0
    if 'MAX_SVD_STEP' in Dictionary:
      val = Dictionary['MAX_SVD_STEP'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # MAX_SVD_STEP must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('1.0\n')

    # INV_ERROR
    check = 0
    if 'INV_ERROR' in Dictionary:
      val = Dictionary['INV_ERROR'][0]
      if val in 'HESSIAN':
        f.write('1\n')
        check = 1
      elif val in 'RF':
        f.write('2\n')
        check = 1
      elif val in 'WORST':
        f.write('3\n')
        check = 1
      else:
        f.write('0\n')
        check = 1
    if check == 0:
      f.write('0\n')

    # PSF_FWHM
    check = 0
    if 'PSF_FWHM' in Dictionary:
      vals = Dictionary['PSF_FWHM']
      # If only one
      if len(vals) == 1:
        # If size 1
        vals[0] = vals[0].split()
        if len(vals[0]) == 1:
          number = False
          try:
            val = interpret(vals[0][0])
            num = float(val)
            number = True
          except ValueError:
            pass
          except:
            verbose(' # PSF_FWHM must be a float, triplets of ' + \
                    'floats, a file, or two floats and a file', \
                    ofolder, verbosity)
            abort(f, filename)
          if number:
            val = interpret(vals[0][0])
            f.write('1\n0\n0\n1d100\n'+val+'\n')
            check = 1
          else:
            f.write('1\n1\n0\n1d100\n'+val+'\n')
            check = 1
      # Not global PSF
      if check == 0:
        NL = len(vals)
        low = []
        up = []
        ff = []
        for val in vals:
          if len(vals) > 1: val = val.split()
          number = False
          try:
            v0 = interpret(val[0])
            v1 = interpret(val[1])
            v2 = interpret(val[2])
            l = float(v0)
            u = float(v1)
            c = float(v2)
            number = True
          except ValueError:
            pass
          except:
            verbose(' # PSF_FWHM must be a float, triplets of ' + \
                    'floats, a file, or two floats and a file', \
                    ofolder, verbosity)
            abort(f, filename)
          if number:
            try:
              val[0] = interpret(val[0])
              val[1] = interpret(val[1])
              val[2] = interpret(val[2])
              l = float(val[0])
              u = float(val[1])
              c = float(val[2])
              if l > u:
                low.append(u)
                up.append(l)
                ff.append(c)
              elif u > l:
                low.append(l)
                up.append(u)
                ff.append(c)
              else:
                verbose(' # PSF_FWHM ranges must have size ' + \
                        'larger than 0 nm', ofolder, verbosity)
                abort(f, filename)
            except:
              verbose(' # PSF_FWHM must be a float or ' + \
                      'triplets of floats', ofolder, verbosity)
              abort(f, filename)
          else:
            try:
              val[0] = interpret(val[0])
              val[1] = interpret(val[1])
              l = float(val[0])
              u = float(val[1])
              c = val[2]
              if l > u:
                low.append(u)
                up.append(l)
                ff.append(c)
              elif u > l:
                low.append(l)
                up.append(u)
                ff.append(c)
              else:
                verbose(' # PSF_FWHM ranges must have size ' + \
                        'larger than 0 nm', ofolder, verbosity)
                abort(f, filename)
            except:
              verbose(' # PSF_FWHM must be a float or ' + \
                      'triplets of floats', ofolder, verbosity)
              abort(f, filename)
        valid,n,low,up,ff = unique_ranges(NL,low,up,ff,[0.,1e100])
        if valid:
          if n > 0:
            f.write('{0}\n'.format(n))
            for l,u,c in zip(low,up,ff):
                if isinstance(c,float):
                    f.write('0\n{0}\n{1}\n{2}\n'.format(l,u,c))
                else:
                    f.write('1\n{0}\n{1}\n{2}\n'.format(l,u,c))
            check = 1
          else:
            f.write('0\n')
            check = 1
        else:
          verbose(' # PSF_FWHM ranges must not intersect', \
                  ofolder, verbosity)
          abort(f, filename)
    if check == 0:
      f.write('0\n')

    # Default bounds
    bounds = {'T': [3500,25000], \
              'VZ': [-20,20], \
              'VT': [0,40], \
              'PG': [0.1,15], \
              'F': [0.,0.95], \
              'J21R': [-1,1], \
              'J21I': [-1,1], \
              'J22R': [-1,1], \
              'J22I': [-1,1]}
    if btyp == 0:
        bounds['B'] = [0,3000]
        bounds['BT'] = [0,math.pi]
        bounds['BP'] = [0,2*math.pi]
    else:
        bounds['B'] = [-3000,3000]
        bounds['BT'] = [0,3000]
        bounds['BP'] = [0,2*math.pi]
    if vtyp == 0:
        bounds['VX'] = [-20,20]
        bounds['VY'] = [-20,20]
    else:
        bounds['VX'] = [0,20]
        bounds['VY'] = [0,2*math.pi]

    # BOUNDS_X
    # For each variable
    for var in varis:
      check = 0
      if 'BOUNDS_'+var in Dictionary:
        val = Dictionary['BOUNDS_'+var]
        try:
          val[0] = interpret(val[0])
          val[1] = interpret(val[1])
          num0 = float(val[0])
          num1 = float(val[1])
          if num1 > num0:
            f.write(val[0] + ' ' + val[1] + '\n')
            check = 1
          elif num1 < num0:
            f.write(val[1] + ' ' + val[0] + '\n')
            check = 1
          else:
            msg = ' # BOUNDS_'+var+' must be '
            msg += 'a list of two unique floats'
            verbose(msg, ofolder, verbosity)
            abort(f, filename)
        except:
          msg = ' # BOUNDS_'+var+' must be '
          msg += 'a list of two unique floats'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        f.write('{0} {1}\n'.format(bounds[var][0], \
                                   bounds[var][1]))

    # EBOUNDS_X
    found = False
    for var in varis:
      if 'EBOUNDS_'+var in Dictionary:
        found = True
        break
    # If found them
    if found:
      out = ['1\n']
      # For each variable
      for var in varis:
        check = 0
        if 'EBOUNDS_'+var in Dictionary:
          vals = Dictionary['EBOUNDS_'+var]
          lout = []
          ncount = 0
          for val in vals:
            try:
              val[0] = interpret(val[0])
              val[1] = interpret(val[1])
              val[2] = interpret(val[2])
              val[3] = interpret(val[3])
              num0 = float(val[0])
              num1 = float(val[1])
              num2 = float(val[2])
              num3 = float(val[3])
              if num1 > num0 and num3 > num2:
                lout.append(val[0] + ' ' + val[1] + ' ' + \
                            val[2] + ' ' + val[3] + '\n')
                ncount += 1
              elif num1 > num0 and num3 < num2:
                lout.append(val[0] + ' ' + val[1] + ' ' + \
                            val[3] + ' ' + val[2] + '\n')
                ncount += 1
              elif num1 < num0 and num3 > num2:
                lout.append(val[1] + ' ' + val[0] + ' ' + \
                            val[2] + ' ' + val[3] + '\n')
                ncount += 1
              elif num1 < num0 and num3 < num2:
                lout.append(val[1] + ' ' + val[0] + ' ' + \
                            val[3] + ' ' + val[2] + '\n')
                ncount += 1
              else:
                msg = ' # EBOUNDS_'+var+' must be '
                msg += 'a list of two pairs of unique floats'
                verbose(msg, ofolder, verbosity)
                abort(f, filename)
            except:
              msg = ' # EBOUNDS_'+var+' must be '
              msg += 'a list of two pairs of unique floats'
              verbose(msg, ofolder, verbosity)
              abort(f, filename)
          if ncount > 0:
            out.append('{0}\n'.format(ncount))
            for llout in lout:
              out.append(llout)
          else:
            out.append('0\n')
          check = 1
        # No boundary for this variable
        if check == 0:
          out.append('0\n')
      # Dump output
      for lout in out:
        f.write(lout)
      del out
    # No special boundaries
    else:
      f.write('0\n')

    # Default scale
    scales = {'B': 100., \
              'BP': 1.8, \
              'T': 2e3, \
              'VZ': 10., \
              'VX': 10., \
              'VT': 10., \
              'PG': 2e0, \
              'F': 1e0, \
              'J21R': 1e-2, \
              'J21I': 1e-2, \
              'J22R': 1e-2, \
              'J22I': 1e-2}
    if btyp == 0:
        scales['BT'] = 1.8
    else:
        scales['BT'] = 100.
    if vtyp == 0:
        scales['VY'] = 10.
    else:
        scales['VY'] = 1.8

    # SCALE_X
    # For each variable
    for var in varis:
      check = 0
      if 'SCALE_'+var in Dictionary:
        val = Dictionary['SCALE_'+var][0]
        try:
          val = interpret(val)
          num = float(val)
          f.write('{0}\n'.format(num))
          check = 1
        except:
          msg = ' # SCALE_'+var+' must be a float'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        f.write('{0}\n'.format(scales[var]))

    # Default perturbation
    pertb = {'B': 1., \
             'BP': 0.03, \
             'T': 50., \
             'VZ': 0.3, \
             'VX': 0.3, \
             'VT': 0.3, \
             'PG': 0.05, \
             'F': 0.01, \
             'J21R': 1e-4, \
             'J21I': 1e-4, \
             'J22R': 1e-4, \
             'J22I': 1e-4}
    if btyp == 0:
        pertb['BT'] = 0.01
    else:
        pertb['BT'] = 2.
    if vtyp == 0:
        pertb['VY'] = 0.3
    else:
        pertb['VY'] = 0.03

    # PERTURB_X
    # For each variable
    for var in varis:
      check = 0
      if 'PERTURB_'+var in Dictionary:
        val = Dictionary['PERTURB_'+var][0]
        try:
          val = interpret(val)
          num = float(val)
          f.write('{0}\n'.format(num))
          check = 1
        except:
          msg = ' # PERTURB_'+var+' must be a float'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        f.write('{0}\n'.format(pertb[var]))

    # MIN_REL_PERTURB_X
    # For each variable
    for var in varis:
      check = 0
      if 'MIN_REL_PERTURB_'+var in Dictionary:
        val = Dictionary['MIN_REL_PERTURB_'+var][0]
        try:
          val = interpret(val)
          num = float(val)
          f.write('{0}\n'.format(num))
          check = 1
        except:
          msg = ' # MIN_REL_PERTURB_'+var+' must be a float'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      if check == 0:
        f.write('0d0\n')

    # INI_BPOS
    check = 0
    if 'INI_BPOS' in Dictionary:
      val = Dictionary['INI_BPOS'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # INI_BPOS must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      if btyp == 0:
        f.write('0.17\n')
      else:
        f.write('10.\n')

    # INI_BAZI
    check = 0
    if 'INI_BAZI' in Dictionary:
      val = Dictionary['INI_BAZI'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # INI_BAZI must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('1.5\n')

    # INI_VPOS
    check = 0
    if 'INI_VPOS' in Dictionary:
      val = Dictionary['INI_VPOS'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # INI_VPOS must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('5.\n')

    # INI_VAZI
    check = 0
    if 'INI_VAZI' in Dictionary:
      val = Dictionary['INI_VAZI'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        verbose(' # INI_VAZI must be a float', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('1.5\n')

    # INV_FRACTION
    check = 0
    if 'INV_FRACTION' in Dictionary:
      val = Dictionary['INV_FRACTION'][0]
      if val in 'YES' or val in 'SI':
        f.write('Y\n')
        check = 1
      elif val in 'NO':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

    # INV_TAU_RANG
    check = 0
    if 'INV_TAU_RANG' in Dictionary:
      val = Dictionary['INV_TAU_RANG']
      try:
        val[0] = interpret(val[0])
        val[1] = interpret(val[1])
        num0 = float(val[0])
        num1 = float(val[1])
        if num1 > num0:
          f.write(val[0] + ' ' + val[1] + '\n')
          check = 1
          tauran = [num0,num1]
        elif num1 < num0:
          f.write(val[1] + ' ' + val[0] + '\n')
          tauran = [num1,num0]
          check = 1
        else:
          msg = ' # INV_TAU_RANG must be '
          msg += 'a list of two unique floats'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      except:
        msg = ' # INV_TAU_RANG must be '
        msg += 'a list of two unique floats'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('-8 1\n')
      tauran = [-8.,1.]

    # ATMO_STRAT
    check = 0
    if 'ATMO_STRAT' in Dictionary:
      vals = Dictionary['ATMO_STRAT']
      NL = len(vals)
      low = []
      up = []
      ff = []
      for val in vals:
        try:
          val[0] = interpret(val[0])
          val[1] = interpret(val[1])
          val[2] = interpret(val[2])
          l = float(val[0])
          u = float(val[1])
          c = float(val[2])
          if c <= 0.:
            verbose(' # ATMO_STRAT third float must be positive', \
                    ofolder, verbosity)
            abort(f, filename)
          if l > u:
            low.append(u)
            up.append(l)
            ff.append(c)
          elif u > l:
            low.append(l)
            up.append(u)
            ff.append(c)
          else:
            verbose(' # ATMO_STRAT ranges must have size ' + \
                    'larger than 0 units', ofolder, verbosity)
            abort(f, filename)
        except:
          verbose(' # ATMO_STRAT wrong format', \
                  ofolder, verbosity)
          abort(f, filename)
      valid,n,low,up,ff = unique_ranges(NL,low,up,ff,tauran)
      if valid:
        if n > 0:
          f.write('{0}\n'.format(n))
          for l,u,c in zip(low,up,ff):
              f.write('{0} {1} {2}\n'.format(l,u,c))
          check = 1
        else:
          f.write('0\n')
          check = 1
      else:
        verbose(' # ATMO_STRAT ranges must not intersect', \
                ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0\n')

    # BROYDEN_LM
    check = 0
    if 'BROYDEN_LM' in Dictionary:
      val = Dictionary['BROYDEN_LM'][0]
      if val in 'YES' or val in 'SI':
        f.write('Y\n')
        check = 1
      elif val in 'NO':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

    # LM_METHOD
    check = 0
    if 'LM_METHOD' in Dictionary:
      val = Dictionary['LM_METHOD'][0]
      if val in 'TRADITIONAL':
        f.write('0\n')
        check = 1
      elif val in 'BACKTRACKING':
        f.write('1\n')
        check = 1
    if check == 0:
      f.write('1\n')

    # LM_BACKTRACKING_MODE
    check = 0
    if 'LM_BACKTRACKING_MODE' in Dictionary:
      val = Dictionary['LM_BACKTRACKING_MODE'][0]
      if val in 'DESPERATION':
        f.write('1\n')
        check = 1
      else:
        f.write('0\n')
        check = 1
    if check == 0:
      f.write('0\n')

    # LM_LAM_BIG_TEST
    check = 0
    if 'LM_LAM_BIG_TEST' in Dictionary:
      val = Dictionary['LM_LAM_BIG_TEST'][0]
      try:
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_BIG_TEST must be a float'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0.1\n')

    # LM_LAM_SMALL_TEST
    check = 0
    if 'LM_LAM_SMALL_TEST' in Dictionary:
      val = Dictionary['LM_LAM_SMALL_TEST'][0]
      try:
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_SMALL_TEST must be a float'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('10.\n')

    # LM_LAM_BIG_PROVE
    check = 0
    if 'LM_LAM_BIG_PROVE' in Dictionary:
      val = Dictionary['LM_LAM_BIG_PROVE'][0]
      try:
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_BIG_PROVE must be a float'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('100.\n')

    # LM_LAM_SMALL_PROVE
    check = 0
    if 'LM_LAM_SMALL_PROVE' in Dictionary:
      val = Dictionary['LM_LAM_SMALL_PROVE'][0]
      try:
        num = float(val)
        f.write(val+'\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_SMALL_PROVE must be a float'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0.1\n')

    # LM_LAMBDA_TRACK
    check = 0
    if 'LM_LAMBDA_TRACK' in Dictionary:
      val = Dictionary['LM_LAMBDA_TRACK'][0]
      try:
        num = int(val)
        f.write(val+'\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_TRACK must be an integer'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('0\n')

    # LM_LAMBDA_RANG
    check = 0
    if 'LM_LAMBDA_RANG' in Dictionary:
      val = Dictionary['LM_LAMBDA_RANG']
      try:
        val[0] = interpret(val[0])
        val[1] = interpret(val[1])
        num0 = float(val[0])
        num1 = float(val[1])
        if num1 > num0:
          f.write(val[0] + ' ' + val[1] + '\n')
          check = 1
        elif num1 < num0:
          f.write(val[1] + ' ' + val[0] + '\n')
          check = 1
        else:
          msg = ' # LM_LAMBDA_RANG must be '
          msg += 'a list of two unique floats'
          verbose(msg, ofolder, verbosity)
          abort(f, filename)
      except:
        msg = ' # LM_LAMBDA_RANG must be '
        msg += 'a list of two unique floats'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('1d-5 1d3\n')

    # LM_LAMBDA_ACCEPT
    check = 0
    if 'LM_LAMBDA_ACCEPT' in Dictionary:
      val = Dictionary['LM_LAMBDA_ACCEPT'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val + '\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_ACCEPT must be a float'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('5d0\n')

    # LM_LAMBDA_REJECT
    check = 0
    if 'LM_LAMBDA_REJECT' in Dictionary:
      val = Dictionary['LM_LAMBDA_REJECT'][0]
      try:
        val = interpret(val)
        num = float(val)
        f.write(val + '\n')
        check = 1
      except:
        msg = ' # LM_LAMBDA_REJECT must be a float'
        verbose(msg, ofolder, verbosity)
        abort(f, filename)
    if check == 0:
      f.write('5d0\n')

    # INV_B_PROJECTION
    check = 0
    if 'INV_B_PROJECTION' in Dictionary:
      val = Dictionary['INV_B_PROJECTION'][0]
      if val in 'YES' or val in 'SI':
        f.write('Y\n')
        check = 1
      elif val in 'NO':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

    # RF_INITSOL
    check = 0
    if 'RF_INITSOL' in Dictionary:
      val = Dictionary['RF_INITSOL'][0]
      if val in 'YES' or val in 'SI':
        f.write('Y\n')
        check = 1
      elif val in 'NO':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('Y\n')

    # INV_NEGL_SIGMA
    check = 0
    if 'INV_NEGL_SIGMA' in Dictionary:
      val = Dictionary['INV_NEGL_SIGMA'][0]
      if val in 'YES' or val in 'SI':
        f.write('Y\n')
        check = 1
      elif val in 'NO':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

    # KEEP_RF
    check = 0
    if 'KEEP_RF' in Dictionary:
      val = Dictionary['KEEP_RF'][0]
      if val in 'YES' or val in 'SI':
        f.write('Y\n')
        check = 1
      elif val in 'NO':
        f.write('N\n')
        check = 1
    if check == 0:
      f.write('N\n')

    # STOREINV_STEP
    check = 0
    valid = 0
    if 'STOREINV_STEP' in Dictionary:
      check = 1
      val = Dictionary['STOREINV_STEP'][0]
      if val.isdigit():
        val = int(val)
        if val > 0:
          f.write('{0:7d}\n'.format(val))
          valid = 1
    if valid == 0 and check == 1:
      f.write('{0:7d}\n'.format(-1))
    if check == 0:
      f.write('{0:7d}\n'.format(-1))

    # FORCE_OBS_FREQ
    check = 0
    if 'FORCE_OBS_FREQ' in Dictionary:
      val = Dictionary['FORCE_OBS_FREQ'][0]
      if val == 'Y' or val == 'YE' or val == 'YES' or \
         val == 'S' or val =='SI':
        f.write('Y\n')
        check = 1
      if val == 'N' or val == 'NO' or val == 'NON':
        f.write('N\n')
        check = 1
      if check == 0:
        f.write('N\n')
    else:
      f.write('N\n')

  f.close()

if __name__ == "__main__":
    rInput()
