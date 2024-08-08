import sys,os,copy,struct
import numpy as np
import scipy.interpolate as scipyinterpolate

debug = False
CS = scipyinterpolate.CubicSpline

######################################################################
######################################################################
######################################################################

def local_help():
    ''' Prints info about this file
    '''

    print('\nThis code creates a Result file combining a full '+ \
          'Result file with a partial newer Result file.')
    print('The use case is to update an old file with potential' + \
          'ongoing inversions that were not finished.')
    print('The combined file allows for using the preliminary ' + \
          'result of incomplete inversions to continue them')
    print('The combined file will never have information on ' + \
          'respose functions')
    print('\n  Usage:')
    print('    python path_to_old_result  ' + \
               'path_to_newer_result  ' + \
               'path_to_combined_result\n')
    sys.exit()

######################################################################
######################################################################
######################################################################

def get_backup(filename):
    ''' Make a copy of the file
    '''

    # Verbose
    print('The file:')
    print('  ',filename)

    # Get new name
    newfilename = filename
    while os.path.isfile(newfilename) or os.path.isdir(newfilename):
        newfilename += '_backup'

    # Verbose
    print('will be overwritten, generating backup:')
    print('  ',newfilename)

    # Copy
    maxbyt = 1000000
    fsize = os.path.getsize(filename)
    fi = open(filename,'rb')
    fo = open(newfilename,'wb')
    while fsize > maxbyt:
        fo.write(fi.read(maxbyt))
        fsize -= maxbyt
    if fsize > 0:
        fo.write(fi.read(fsize))
    fi.close()
    fo.close()

    return newfilename

######################################################################
######################################################################
######################################################################

def are_different(x,y,e):
    ''' Check if two arrays are different
    '''
    if len(x.shape) != len(y.shape): return True
    for m,n in zip(x.shape,y.shape):
        if m != n: return True
    if np.max(np.absolute(x.flatten() - y.flatten())/ \
              np.absolute(x.flatten() + y.flatten())) > 0.5*e: \
       return True
    return False

######################################################################
######################################################################
######################################################################

def are_compatible(h1,h2):
    ''' Check if the two files to combine are compatible
    '''
    if h1['nx'] != h2['nx']:
        print('Both files have different X dimensions and thus ' + \
                'are not compatible:',h1['nx'],'!=',h2['nx'])
        return False
    if h1['ny'] != h2['ny']:
        print('Both files have different Y dimensions and thus ' + \
                'are not compatible:',h1['ny'],'!=',h2['ny'])
        return False
    if h1['nl'] != h2['nl']:
        print('Both files have different wavelength dimensions ' + \
              'and thus are not compatible:',h1['nl'],'!=',h2['nl'])
        return False
    if h1['nvar'] != h2['nvar']:
        print('Both files have different number of inverted ' + \
              'variables and thus are not compatible:', \
              h1['nvar'],'!=',h2['nvar'])
        return False
    if h1['polarization'] != h2['polarization']:
        if h1['jkqa']:
            print('The original file has polarization data, but ' + \
                  'the newer one does not. The combined file ' + \
                  'will drop the polarization')
        else:
            print('The newer file has polarization data, but the' + \
                  'original newer one does not. The combined ' + \
                  'file will have zero polarization wher there ' + \
                  'is no newer data')
    if h1['jkqa'] != h2['jkqa']:
        if h1['jkqa']:
            print('The original file has JKQ data, but the newer ' + \
                  'one does not. The combined file will fill ' + \
                  'with zeros the columns with newer results')
        else:
            print('The newer file has JKQ data, but the original ' + \
                  'one does not. The combined file will fill ' + \
                  'with zeros the columns without newer results')
    if h1['vtype'] != h2['vtype']:
        print('The files have different type of velocity data. ' + \
              'Because there is no LOS information in these ' + \
              'files, the combined file will be in vertical')
    if h1['btype'] != h2['btype']:
        print('The files have different type of magnetic data. ' + \
              'Because there is no LOS information in these ' + \
              'files, the combined file will be in vertical')
    return True


######################################################################
######################################################################
######################################################################

def get_head(filename):
    ''' Extract the header data
    '''

    # Initialize
    head = {'check': False}

    # Open
    f = open(filename,'rb')

    # Label
    label = f.read(4).decode('utf-8')
    if label != 'invo':
        print('Wrong file label. Expected "invo", got: ',label)
        f.close()
        return head

    # Read metadata inversion
    head['info'] = struct.unpack('i',f.read(4))[0]

    # Interpretation
    info = head['info']
    head['jkqa'] = info > 7
    if info > 7: info -= 8
    if info > 3:
        head['vtype'] = 1
        info -= 4
    else:
        head['vtype'] = 0
    # B type
    if info > 1:
        head['btype'] = 1
        info -= 2
    else:
        head['btype'] = 0
    # Polarization
    head['polarization'] = info > 0
    # Dimensions
    head['nx'] = struct.unpack('i',f.read(4))[0]
    head['ny'] = struct.unpack('i',f.read(4))[0]
    head['nz'] = struct.unpack('i',f.read(4))[0]

    #
    # Prepare jumps atmosphere related

    # Header
    head['s_head'] = 5*4

    # Jump to beginning of atmosphere
    head['jump_to_atmo'] = head['s_head']

    # Number of variables in atmospheric model
    if head['jkqa']:
      head['nvar_atmo'] = 27
    else:
      head['nvar_atmo'] = 19

    # Atmosphere pixel size                             v f_diff
    head['s_atmo_c'] = head['nz']*head['nvar_atmo']*4 + 4

    # Atmosphere size
    head['s_atmo'] = head['nx']*head['ny']*head['s_atmo_c']

    # Skip atmosphere
    f.seek(head['s_atmo'],1)

    # Get header for result
    head['nvar'] = struct.unpack('i',f.read(4))[0]
    head['nl'] = struct.unpack('i',f.read(4))[0]
    head['jump_to_lambda'] = head['jump_to_atmo'] + head['s_atmo'] + 8
    head['lambda'] = np.array(struct.unpack('d'*head['nl'], \
                                            f.read(8*head['nl'])))

    # Prepare lists with flags and nodes for inversion
    # variables
    head['inv_flag'] = []
    head['inv_node'] = []

    # For each variable
    for ivar in range(head['nvar']):

        # Read if inverting
        head['inv_flag'].append(struct.unpack('i', f.read(4))[0] > 0)

        # Number of nodes
        head['inv_node'].append(struct.unpack('i', f.read(4))[0])

    #
    # Prepare jumps results related

    # Size of header
    head['s_res_h'] = 4 + 4 + head['nvar']*8 + head['nl']*8

    # Size of results, Stokes part
    if head['polarization']:

        head['s_res_c'] = 8 + 32*head['nl']

    else:

        head['s_res_c'] = 8 + 8*head['nl']

    # For each variable, add nodes info
    for flag,node in zip(head['inv_flag'],head['inv_node']):

        # Nodes
        if node <= 0: continue

        # Inverting
        if flag:
            head['s_res_c'] += 12*node
        else:
            head['s_res_c'] += 8*node

    # Results size
    head['s_res'] = head['nx']*head['ny']*head['s_res_c']

    # Jump to beginning of results
    head['jump_to_res'] = head['s_head'] + head['s_atmo'] + \
                          head['s_res_h']

    # Validate and return
    head['check'] = True
    head['f'] = filename

    # Close and return
    f.close()
    return head

######################################################################
######################################################################
######################################################################

def init_file(h):
    ''' Create and fill header for the new combined file
    '''

    # Open destiny file
    try:
        f = open(h['f'],'wb')
    except:
        raise

    # Label
    f.write(h['label'].encode())

    # Info
    f.write(struct.pack('i',h['info']))

    # Dimensions
    f.write(struct.pack('i',h['nx']))
    f.write(struct.pack('i',h['ny']))
    f.write(struct.pack('i',h['nz']))

    # Header
    h['s_head'] = 5*4

    # Jump to beginning of atmosphere
    h['jump_to_atmo'] = h['s_head']

    # Number of variables in atmospheric model
    if h['jkqa']:
      h['nvar_atmo'] = 27
    else:
      h['nvar_atmo'] = 19

    # Atmosphere pixel size                    v f_diff
    h['s_atmo_c'] = h['nz']*h['nvar_atmo']*4 + 4

    # Atmosphere size
    h['s_atmo'] = h['nx']*h['ny']*h['s_atmo_c']


    # Close and return
    f.close()
    return h

######################################################################
######################################################################
######################################################################

def step2_file(h):
    ''' Continue writing file between atmo and results
    '''

    # Open destiny file
    try:
        f = open(h['f'],'ab')
    except:
        raise

    # Skip atmosphere
    f.seek(h['jump_to_atmo'] + h['s_atmo'],0)

    # Get header for result
    f.write(struct.pack('i',h['nvar']))
    f.write(struct.pack('i',h['nl']))
    h['jump_to_lambda'] = h['jump_to_atmo'] + h['s_atmo'] + 8
    f.write(struct.pack('d'*h['nl'],*h['lambda']))

    # For each variable
    for ivar in range(h['nvar']):

        # Read if inverting
        if h['inv_flag'][ivar]:
            f.write(struct.pack('i',1))
        else:
            f.write(struct.pack('i',0))

        # Number of nodes
        f.write(struct.pack('i',h['inv_node'][ivar]))

    #
    # Prepare jumps results related

    # Size of header
    h['s_res_h'] = 4 + 4 + h['nvar']*8 + h['nl']*8

    # Size of results, Stokes part
    if h['polarization']:

        h['s_res_c'] = 8 + 32*h['nl']

    else:

        h['s_res_c'] = 8 + 8*h['nl']

    # For each variable, add nodes info
    for flag,node in zip(h['inv_flag'],h['inv_node']):

        # Nodes
        if node <= 0: continue

        # Inverting
        if flag:
            h['s_res_c'] += 12*node
        else:
            h['s_res_c'] += 8*node

    # Results size
    h['s_res'] = h['nx']*h['ny']*h['s_res_c']

    # Jump to beginning of results
    h['jump_to_res'] = h['s_head'] + h['s_atmo'] + h['s_res_h']

    # Close and return
    f.close()
    return h

######################################################################
######################################################################
######################################################################

def set_atmos(h1,h2,h3):
    ''' Fill atmosphere
    '''

    # Open files and seek to atmo
    f1 = open(h1['f'],'rb')
    f1.seek(h1['jump_to_atmo'],0)
    f2 = open(h2['f'],'rb')
    f2.seek(h2['jump_to_atmo'],0)
    f3 = open(h3['f'],'ab')
    f3.seek(h3['jump_to_atmo'],0)

    # Sizes
    bsiz1 = h1['s_atmo_c']
    siz1 = bsiz1//4
    bsiz2 = h2['s_atmo_c']
    siz2 = bsiz2//4
    siz3 = h3['nz']*h3['nvar']

    # Newer
    newer = np.zeros((h3['nx'],h3['ny']))

    # For each column
    for ix in range(h1['nx']):
        for iy in range(h1['ny']):

            # Get data from file origin
            d1 = np.array(struct.unpack('f'*(siz1-1), \
                                        f1.read(bsiz1-4))). \
                         reshape((h1['nvar_atmo'],h1['nz']))
            fd1 = struct.unpack('f',f1.read(4))[0]

            # Get data from file newer
            d2 = np.array(struct.unpack('f'*(siz2-1), \
                                        f2.read(bsiz2-4))). \
                         reshape((h2['nvar_atmo'],h2['nz']))
            fd2 = struct.unpack('f',f2.read(4))[0]

            # Valid newer?
            valid = True

            # Check positive
            if np.min(d1[0,:]) <= 0.:
                valid = False

            # Check monotonic
            if np.min(d1[0,1:] - d1[0,:-1]) <= 0.:
                valid = False

            # Check reasonable
            if np.max(d1[0,:]) < 1e-10 or np.min(d1[0,:]) > 1e10:
                valid = False

            # Newer data
            if valid:

                # Set new
                newer[ix,iy] = 1

                # If different
                if h3['ldiff']:

                    # Log tau
                    l1 = np.log10(d1[0,:])
                    l2 = np.log10(d2[0,:])

                    # Create new array
                    d3 = np.empty((h3['nvar_atmo'],h3['nz']), \
                                  dtype=np.float32)

                    # Number of elements
                    siz3 = h3['nvar_atmo']*h3['nz']

                    # Copy heights
                    d3[0,:] = d1['tau']

                    # For all variables not JKQ
                    for i in range(1,19):

                        # If constant
                        if np.min(d2[i,:]) >= np.max(d2[i,:]):

                            # Copy just one value
                            d3[0,:] = d2[i,0]

                        # Not constant
                        else:

                            # Interpolate
                            d3[i,:] = CS(l2,d2[i,:])(l1)

                    # If there are JKQ
                    if h3['jkqa']:

                        # If there are in newer
                        if h2['jkqa']:

                            # For all variables not JKQ
                            for i in range(20,27):

                                # If constant
                                if np.min(d2[i,:]) >= np.max(d2[i,:]):

                                    # Copy just one value
                                    d3[0,:] = d2[i,0]

                                # Not constant
                                else:

                                    # Interpolate
                                    d3[i,:] = CS(l2,d2[i,:])(l1)

                        # If there are not in newer
                        else:

                            d3[20:27,:] = np.zeros((8,h3['nz']), \
                                                   dtype=np.float32)

                    # Write new data
                    f3.write(struct.pack('f'*(siz3-1), \
                                         *(d3.flatten())))

                # Same axis
                else:

                    # Write new data
                    f3.write(struct.pack('f'*(siz2-1), \
                                         *(d2.flatten())))
                    # If fill with JKQ but there is not
                    if h3['jkqa'] and not h2['jkqa']:
                        f3.write(struct.pack('f'*h2['nz']*8, \
                                             np.zeros((h2['nz']*8), \
                                                      dtype=np.float32)))
                # Diffuse light
                f3.write(struct.pack('f',fd2))

            # Old data
            else:

                # Write old data
                f3.write(struct.pack('f'*(siz1-1), \
                                     *(d1.flatten())))
                # If fill with JKQ but there is not
                if h3['jkqa'] and not h1['jkqa']:
                    f3.write(struct.pack('f'*h1['nz']*8, \
                                         np.zeros((h1['nz']*8), \
                                         dtype=np.float32)))
                # Diffuse light
                f3.write(struct.pack('f',fd1))

    # Close files
    f1.close()
    f2.close()
    f3.close()

    # Return
    return newer

######################################################################
######################################################################
######################################################################

def set_nodes(h1,h2,h3,newer):
    ''' Fill atmosphere
    '''

    # Get column
    def get_column_atmo(h,ix,iy):
        o = {}
        f = open(h['f'],'rb')
        bsiz = h['s_atmo_c']
        siz = bsiz//4
        f = open(h['f'],'rb')
        f.seek(h['jump_to_atmo'] + iy*bsiz + h['ny']*ix*bsiz,0)
        bsiz -= 4
        siz -= 1
        o['data'] = np.array(struct.unpack('f'*siz, \
                                         f.read(bsiz))). \
                           reshape((h['nvar_atmo'],h['nz']))
        o['f_diff'] = struct.unpack('f',f.read(4))[0]
        f.close()
        return o

    # Get column
    def get_column(h,ix,iy):
        o = {}
        f = open(h['f'],'rb')
        bsiz = h['s_res_c']
        siz = bsiz//4
        f.seek(h['jump_to_res'] + iy*bsiz + h['ny']*ix*bsiz,0)
        o['chi0'] = struct.unpack('f',f.read(4))[0]
        o['chi'] = struct.unpack('f',f.read(4))[0]
        o['stko'] = np.zeros((4,h['nl']),dtype=np.float32)
        o['stkf'] = np.zeros((4,h['nl']),dtype=np.float32)
        if h['polarization']:
            o['stko'] = np.array(struct.unpack('f'*h['nl']*4, \
                                               f.read(16*h['nl']))). \
                                reshape((4,h['nl']))
            o['stkf'] = np.array(struct.unpack('f'*h['nl']*4, \
                                               f.read(16*h['nl']))). \
                                reshape((4,h['nl']))
        else:
            o['stko'][0,:] = np.array(struct.unpack('f'*h['nl'], \
                                               f.read(4*h['nl'])))
            o['stkf'][0,:] = np.array(struct.unpack('f'*h['nl'], \
                                               f.read(4*h['nl'])))
        for ipar,node in enumerate(h['inv_node']):
            if node <= 0: continue
            o[ipar] = np.empty((3,node),dtype=np.float32)
            o[ipar][0:2,:] = np.array(struct.unpack('f'*node*2, \
                                                    f.read(8*node))). \
                                      reshape((2,node))
            if h['inv_flag'][ipar]:
                o[ipar][2,:] = np.array(struct.unpack('f'*node, \
                                                      f.read(4*node)))
            else:
                o[ipar][2,:] = 0e0
        f.close()
        return o

    # Open for nodes
    f3 = open(h3['f'],'ab')
    f3.seek(h3['jump_to_res'],0)

    # Dictionary of variables
    kvar = {0:  3, \
            1:  4, \
            2:  5, \
            3: -1, \
            4:  2, \
            5:  6, \
            6:  7, \
            7:  8, \
            8:  9, \
            9:  2, \
           10: 23, \
           11: 24, \
           12: 25, \
           13: 26}

    # For each column
    for ix in range(h3['nx']):
        for iy in range(h3['ny']):

            #
            # Stokes parameters
            #

            # If new
            if newer[ix,iy]:
                h = h2
            else:
                h = h1

            # Get node data
            node = get_column(h,ix,iy)

            # Write chi
            f3.write(struct.pack('f',node['chi0']))
            f3.write(struct.pack('f',node['chi']))

            # If polarization
            if h3['polarization']:

                # Write
                f3.write(struct.pack('f'*h3['nl']*4, \
                                     *(node['stko'].flatten())))
                f3.write(struct.pack('f'*h3['nl']*4, \
                                     *(node['stkf'].flatten())))
            # No polarization
            else:

                # Write
                f3.write(struct.pack('f'*h3['nl'], \
                                     *(node['stko'][0,:])))
                f3.write(struct.pack('f'*h3['nl'], \
                                     *(node['stkf'][0,:])))

            #
            # Nodes
            #

            # Column
            col = get_column_atmo(h,ix,iy)

            # For each variable
            for ivar,n in enumerate(h3['inv_node']):

                # If number of nodes is zero
                if n <= 0: continue

                # If different number of nodes
                if h['inv_node'][ivar] != n:

                    # Get atmos variable
                    jvar = kvar[ivar]

                    # Get space
                    onode = np.zeros((3,n),dtype=np.float32)

                    # Log tau
                    ltau = np.log10(col['data'][0,:])

                    # If pressure or f
                    if jvar == 2 or jvar < 0:
                        onode[ivar][0,0] = ltau[0]
                        if jvar == 2:
                            onode[1,0] = col['data'][jvar,0]
                        else:
                            onode[1,0] = col['f_diff']

                    # Stratified
                    else:

                        # Choose equidistant taus
                        onode[0,:] = np.linspace(ltau[0], \
                                                 ltau[-1], \
                                                 n,endpoint=True)

                    # Stratifications
                    if jvar <= 18 or \
                        (jvar > 18 and col['data'].shape[0] > 19):

                        # Copy data
                        for i in range(n):
                            j = np.argmin( \
                                    np.absolute(ltau - onode[0,i]))
                            onode[0,i] = ltau[j]
                            onode[1,i] = col['data'][jvar,j]

                # Equal
                else:

                    # Copy
                    onode = node[ivar]

                # Write
                f3.write(struct.pack('f'*n,*onode[0,:]))
                f3.write(struct.pack('f'*n,*onode[1,:]))
                if h3['inv_flag'][ivar]:
                    f3.write(struct.pack('f'*n,*onode[2,:]))

    # Close
    f3.close()

######################################################################
######################################################################
######################################################################

def setup_head(h1,h2):
    ''' Prepare the header data for the combined file
    '''

    # Initialize
    h3 = {}

    # Label
    h3['label'] = 'invo'

    # JKQa
    h3['jkqa'] = h1['jkqa'] or h2['jkqa']

    # Type from the original
    for typ in ['vtype','btype']:
        if h1[typ] == h2[typ]:
            h3[typ] = h1[typ]
        else:
            h3[typ] = 0

    # Polarization
    h3['polarization'] = h1['polarization'] or h2['polarization']

    # Dimensions
    h3['nx'] = h1['nx']
    h3['ny'] = h1['ny']
    h3['nz'] = h1['nz']

    # nvar_atmo
    if h3['jkqa']:
        h3['nvar_atmo'] = 27
    else:
        h3['nvar_atmo'] = 19

    # Prepare list of flags and nodes for inversion variables
    h3['nvar'] = h1['nvar']
    h3['nl'] = h1['nl']
    h3['inv_flag'] = []
    h3['inv_node'] = []

    # Run over variables
    for ivar in range(h1['nvar']):
        h3['inv_flag'].append(h1['inv_flag'][ivar] or \
                              h2['inv_flag'][ivar])
        if h3['inv_flag'][ivar]:
            if h1['inv_flag'][ivar]:
                h3['inv_node'].append(h1['inv_node'][ivar])
            else:
                h3['inv_node'].append(h2['inv_node'][ivar])
        else:
            h3['inv_node'].append(0)

    # Information
    h3['info'] = 0
    if h3['polarization']: h3['info'] += 1
    if h3['vtype']: h3['info'] += 2
    if h3['btype']: h3['info'] += 4
    if h3['jkqa']: h3['info'] += 8

    # Copy lambda
    h3['lambda'] = h1['lambda']

    # Return head
    return h3

######################################################################
######################################################################
######################################################################

def get_ltau(h):
    ''' Get the log tau scale
    '''

    # Open file
    f = open(h['f'],'rb')

    # Jump to atmospheric data
    f.seek(h['jump_to_atmo'],0)

    # Initialize
    icol = 0
    ncol = h['nx']*h['ny']

    # Sizes
    csiz = h['nz']
    cbsiz = csiz*4
    bsiz = h['s_atmo_c'] - cbsiz
    siz = bsiz//4

    # For each column
    while True:

        # Get tau
        z = np.array(struct.unpack('f'*csiz, \
                                   f.read(cbsiz)))

        # Initialize valid
        valid = True

        # Check positive
        if np.min(z) <= 0.:
            valid = False

        # Check monotonic
        if np.min(z[1:] - z[:-1]) <= 0.:
            valid = False

        # Check reasonable
        if np.max(z) < 1e-10 or np.min(z) > 1e10:
            valid = False

        # If valid, leave
        if valid: break

        # If not valid, count
        icol += 1

        # Reached limit
        if icol >= ncol:
            f.close()
            print('Cound not find any valid tau scale in file:')
            print('  ',h['f'])
            print('Are you sure there is any data in there?')
            sys.exit()

        # Skip rest of column
        f.seek(bsiz,1)

    # Close
    f.close()

    if debug: print(  'ltau: ',np.log10(z))

    # Return
    return np.log10(z)

######################################################################
######################################################################
######################################################################

def run_combination(h1,h2,filename3):
    ''' Create the new combined file
    '''

    # Setup new header
    h3 = setup_head(h1,h2)
    h3['f'] = filename3

    if debug: print('Header 3:',h3)

    # Get height scales from sources
    ltau1 = get_ltau(h1)
    ltau2 = get_ltau(h2)

    # Are they different?
    ldiff = are_different(ltau1,ltau2,1e-16)

    if debug: print('Different ltau: ',ldiff)

    # Initialize header
    h3 = init_file(h3)
    h3['ldiff'] = ldiff

    # Update atmosphere
    newer = set_atmos(h1,h2,h3)

    # Continue writing
    h3 = step2_file(h3)

    # Update nodes
    set_nodes(h1,h2,h3,newer)

######################################################################
######################################################################
######################################################################

def main():
    ''' Main part of the code
    '''

    # Help needed?
    if len(sys.argv) < 4:
        local_help()
    if sys.argv[1].lower().strip() == 'help':
        local_help()

    # Check path 1
    file1 = sys.argv[1]
    try:
        f = open(file1,'r')
        f.close()
    except FileNotFoundError:
        print('Aborting because could not find the old result file:')
        print('  ',file1)
        sys.exit()
    except:
        raise

    # Check path 2
    file2 = sys.argv[2]
    try:
        f = open(file2,'r')
        f.close()
    except FileNotFoundError:
        print('Aborting because could not find the newer result ' + \
              'file:')
        print('  ',file2)
        sys.exit()
    except:
        raise

    # Check path 3
    file3 = sys.argv[3]
    if file1 == file2:
        print('Aborting because the old and new files are ' + \
              'the same:')
        print('  ',file1)
        sys.exit()
    if file3 == file1: file1 = get_backup(file1)
    if file3 == file2: file2 = get_backup(file2)

    # Get header
    head1 = get_head(file1)
    if not head1['check']: sys.exit()
    head2 = get_head(file2)
    if not head1['check']: sys.exit()

    if debug:
        print('Header 1:',head1)
        print('Header 2:',head2)

    # Compatible?
    if not are_compatible(head1,head2):
        sys.exit()

    # Combine
    run_combination(head1,head2,file3)

    # Verbose
    print('Created file: ')
    print('  ',file3)


######################################################################
######################################################################
######################################################################

if __name__ == '__main__':
    main()
