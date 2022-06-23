import argparse
parser = argparse.ArgumentParser(description= "add N val at beg and end of bed rec")
parser.add_argument("-bed", "--b", help="output Bed file",type=str)
parser.add_argument("-number", "--N", help="N val to be add")
parser.add_argument("-out", "--o", help="output file",type=str)
arg = parser.parse_args()
filehandle = open(arg.b)
oWriter=open(arg.o, "w")
N=int(arg.N)
for line in filehandle:
    stripped=line.split("\t")
    if stripped[3] > stripped[4]:
        end = str(int(stripped[3]) + N)
        start = str(int(stripped[4]) - N)
    else:
        end = str(int(stripped[4]) + N)
        start = str(int(stripped[3]) - N)
    oWriter.write(stripped[0]+ "\t" +start + "\t"+ end + "\t"+ stripped[8])