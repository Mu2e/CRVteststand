#! /usr/bin/python3

# automatic transmission script. Parallelly processing multiple rsyncs at the same time

import glob, os, sys
import time
import subprocess
import multiprocessing as mp

LIMIT_BYTES = 20 * 1024**3 # max file size allowed

gpvmPool = ["mu2egpvm04.fnal.gov", "mu2egpvm05.fnal.gov"]
gpvmUsgCount = mp.Manager().Array('i',[0, 0])
lock = mp.Manager().Lock()
userName = "mu2epro"
#userName = "yongyiwu"
nMaxProcessPerGpvm = min(3, int(mp.cpu_count()/len(gpvmPool)))
#nMaxProcessPerGpvm = min(2, int(mp.cpu_count()/len(gpvmPool)))
destination = "/pnfs/mu2e/persistent/users/mu2epro/RDM/remote/crv/output/"
#destination = "/mu2e/data/crv/wideband/"
print ("CPU count = ", mp.cpu_count())
print ("Destination GPVMs are ", gpvmPool)
print ("Destination dir is ", destination)
print ("Max transmission per GPVM is", nMaxProcessPerGpvm)

fileSpec = sys.argv[1]
filelist = glob.glob(fileSpec)
print (len(filelist), "files to transmit...\n")

for tfile in filelist:
    try:
        size = os.path.getsize(tfile)
    except OSError as e:
        print(f"ERROR: Cannot access '{tfile}': {e}", file=sys.stderr)
        sys.exit(1)
    if size > LIMIT_BYTES:
        print(f"ERROR: File '{tfile}' is too large: {size / (1024**3):.2f} GiB", file=sys.stderr)
        sys.exit(1)

if len(sys.argv) > 2:
    userName = sys.argv[2]

def transmitOneFile(filename):
    global gpvmPool, gpvmUsgCount, userName, destination, filelist
    
    lock.acquire()
    iGpvm = 0
    # too bad there's no numpy and cannot be installed...
    # essentially doing iGpvm = np.argmin(gpvmUsgCount)
    mincount = min(gpvmUsgCount)
    for i in range(len(gpvmUsgCount)):
        if gpvmUsgCount[i] == mincount:
            iGpvm = i
            break      
    tGpvm = gpvmPool[iGpvm]
    print ("============================================================")
    print ("Starting", filename.split('/')[-1], "on", tGpvm, "...")
    print ("GPVM occupied:", gpvmUsgCount, end = ' => ')
    gpvmUsgCount[iGpvm] += 1
    print (gpvmUsgCount)
    sys.stdout.flush()
    lock.release()
    
    t0 = time.time()
    # cmd = "rsync -avzh --progress %s %s@%s:%s"%(filename, userName, tGpvm, destination)
    cmd = "rsync -avzh %s %s@%s:%s"%(filename, userName, tGpvm, destination)
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    lock.acquire()
    print ("============================================================")
    print (filename.split('/')[-1], "transmitting through", tGpvm)
    print (output.decode('ascii'))
    t1 = time.time()
    elapsed = t1-t0
    print ("Time elapsed " + time.strftime("%H:%M:%S.{}".format(str(elapsed % 1)[2:])[:10], time.gmtime(elapsed)))
    print ("Releasing link to", tGpvm)
    print ("GPVM occupied:", gpvmUsgCount, end = ' => ')
    gpvmUsgCount[iGpvm] -= 1
    print (gpvmUsgCount)
    sys.stdout.flush()
    lock.release()
    
    return

def main():
    global gpvmPool, nMaxProcessPerGpvm
    nWorker = nMaxProcessPerGpvm*len(gpvmPool)
    p = mp.Pool(processes=nWorker)
    for file in filelist:
        p.apply_async(transmitOneFile, (file,))
    p.close()
    p.join()
    return

if __name__ == "__main__":
    main()
