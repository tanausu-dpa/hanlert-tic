import struct,sys

def main():

    try:
        f = open(sys.argv[1],'rb')
    except:
        raise
   #label = str(f.read(5)).decode('utf-8')
    label = f.read(5).decode('utf-8')
    print(label)
    nx = int(struct.unpack('i',f.read(4))[0])
    ny = int(struct.unpack('i',f.read(4))[0])
    nn = nx*ny
    ni = 0
    fi = 0
    iv = 0
    while True:
        try:
            c = int(struct.unpack('i',f.read(4))[0])
            if c > 0:
                iv += 1
            elif c < 0:
                fi += 1
            ni += 1
            if ni == nn:
               break
        except:
            break
    print(iv,' of ',nn,'  (',ni,')','  Failed:',fi,'. Total:',iv+fi)
    f.close()

    print('Failed:')
    try:
        f = open(sys.argv[1],'rb')
    except:
        raise
   #label = str(f.read(5)).decode('utf-8')
    label = f.read(5).decode('utf-8')
    nx = int(struct.unpack('i',f.read(4))[0])
    ny = int(struct.unpack('i',f.read(4))[0])
    nn = nx*ny
    ii = 0
    for ix in range(nx):
        for iy in range(ny):
            try:
                c = int(struct.unpack('i',f.read(4))[0])
                ii += 1
                if c < 0:
                    print('x',ix+1,' y',iy+1, 'los',ii)
            except:
                break
    f.close()


if __name__ == '__main__':
    main()
