#! /usr/bin/python3

import sys, io, re, os

IN_BUF_SIZE = 1024*1024*256 # 256 MB
OUT_BUF_SIZE = 1024*1024*512 # 512 MB dumps
d = 1024 # rough # of bytes slightly large than that between 'begin of spill' and 'FEB0'

def read_chunk(file,chunk_size = IN_BUF_SIZE):
    while True:
        data = file.read(chunk_size)
        if not data:
            break
        yield data

infs = sys.argv[1]
outfs = sys.argv[2]
outfs0 = outfs
outfn = outfs.split("/")[-1]
outdir = "/".join(outfs.split("/")[:-1])

nev = int(sys.argv[3])
nsplit = int(sys.argv[4])

# n file expected
nftot = int(nev/nsplit)
# ev per file, remainder goes in last file
nevpf = int(nev/nftot)

# skip this many files/ev (first actual output will have this suffix)
if len(sys.argv)>5:
    nskipf = int(sys.argv[5])
    nskipev = nskipf*nevpf
else:
    nskipev = -1

seq = 0

infile = io.open(infs, "rb", buffering=IN_BUF_SIZE)

nev = 0
nf = 0

beginFileFlag = True
last_data = None

outfile = None
if nskipev < 0:
    outfile = io.open(outfs, "wb", buffering=OUT_BUF_SIZE)
    # outfile.reconfigure(write_through=True) # no need already binary; only apply for text files 
    print("Opening ",outfn)
else:
    outfile = io.open("trash.txt","wb")

for chunk in read_chunk(infile):
    
    data = None
    if last_data:
        data = last_data + chunk
        last_data = None
    else:
        data = chunk
    
    beginEvent = [m.start() for m in re.finditer(b'>\x14--\*\* SOURCE = FEB0\r\n', data)]
    if not beginEvent: # contents all belong to previous spill (event) or header
        if nev > nskipev:
            outfile.write(data[:(-d)]) 
        last_data = data[(-d):]
        # begin of spill and start of FEB0 always appear within 1kB.
        # this treatment prevents start of spill appear in the chunk but start of FEB0 in the next
    else:
        beginSpill = [m.start() for m in re.finditer(b'\x12--Begin of spill\r\n', data)]
        indexSpill = 0
        splitPt = []
        for iEvent, bEvent in enumerate(beginEvent):
            for iSpill, bSpill in enumerate(beginSpill[indexSpill:]):
                if bSpill<bEvent and bSpill>(beginEvent[iEvent-1] if iEvent>0 else -1):
                    if iSpill == 0 or bSpill-tSplitPt>100: # empty spill belongs to the next spill
                        tSplitPt = bSpill
                else:
                    splitPt.append(tSplitPt)
                    indexSpill = beginSpill.index(bSpill)
        if len(splitPt)<len(beginEvent):
            splitPt.append(tSplitPt) # take care of the last spill 
        # at this point, data[:(splitPt[0])] belongs to previous spill / run header
        # data[(splitPt[0]):(splitPt[-1])] goes to the spills
        # if splitPt[-1] is before index (-2d), data[(splitPt[-1]):(-d)] always in a spill and use [(-d):] for stub
        # if splitPt[-1] is or after index (-2d), data[(splitPt[-1]):] goes to stub
        # such scheme ensures there's always a 'begin of spill' before the 'FEB0'
        if data[:(splitPt[0])] and nev > nskipev:
            outfile.write(data[:(splitPt[0])])
        
        if splitPt[-1] < (len(data)-2*d):
            splitPt.append((-d))

        if data[(splitPt[0]):(splitPt[-1])]:
            for i in range(len(splitPt)-1):
                # this spill(event) is data[splitPt[i]:splitPt[i+1]]

                if nev>0 and nev%nevpf==0 and nf!=(nftot-1):
                    nf = nf + 1
                    seq = seq + 1
                    if nev >= nskipev: # close file and start new
                        outfile.close()
                        sseq = "_{:03d}.".format(seq)
                        outfs = outfs0.replace("_000.",sseq)
                        outfn = outfs.split("/")[-1]
                        outfile = io.open(outfs, "wb", buffering=OUT_BUF_SIZE)
                        # outfile.reconfigure(write_through=True) # no need already binary; only apply for text files 
                        print ("Opening ", outfn)

                nev = nev + 1
                if nev > nskipev:
                    if nev%200 == 0:
                        print ("Processing event ", nev)
                else:
                    if nev%nsplit == 0:
                        print ("Skipping event ", nev)

                if nev > nskipev:
                    outfile.write(data[splitPt[i]:splitPt[i+1]])

        last_data = data[(splitPt[-1]):]

# take care of whatever is left in last_data 
nExtraEvent = len([m.start() for m in re.finditer(b'>\x14--\*\* SOURCE = FEB0\r\n', last_data)])
nev = nev + nExtraEvent
outfile.write(last_data)
outfile.close()
infile.close()

nf = nf + 1

print ("Files: wrote ", nf)
print ("Events: processed ", nev)

# print ("=== Sanity Check ===")
# cmd1 = "wc -l " + infs
# cmd2 = "wc -l " + outfs0.replace("_000.","_*.")
# print ("Calling: "+cmd1)
# os.system(cmd1)
# print ("Calling: "+cmd2)
# os.system(cmd2)
