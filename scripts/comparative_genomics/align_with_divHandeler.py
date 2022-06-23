import argparse
import math
#Def
def variance(data):
# Number of observations
    n = len(data)
    # Mean of the data
    mean = sum(data) / n
    # Square deviations
    deviations = [(x - mean) ** 2 for x in data]
    # Variance
    variance = sum(deviations) / n
    return variance
def stdev(data):
    var = variance(data)
    std_dev = math.sqrt(var)
    return std_dev



parser = argparse.ArgumentParser(description= "Select Kimura under param and create bed")
parser.add_argument("--awd", "-align_with_div",  help="out of calcDivergenceFromAlign.pl")
parser.add_argument("--k", "-kimura",  help="Kimura threshold")
parser.add_argument("-out", "--o", help="output file",type=str)
parser.add_argument("-bed", "--b", help="output Bed file",type=str)
arg = parser.parse_args()
filehandle = open(arg.awd)
Tranver={}
name={}
Id="asd"
transitions_list=[]
for line in filehandle:
        if line.startswith(" "):
                continue
        elif line.startswith("\t"):
                continue
        elif line.startswith("Matrix"):
                continue
        elif line.startswith("Gap_init"):
                continue
        elif line.startswith("C"):
                continue
        elif line[0].isdigit():
            stripped=line.split(" ")
            if  stripped[8].startswith("C"):
                TE=stripped[9]
            else:
                TE=stripped[8]
            Id = (stripped[4]+"\t"+ stripped[5]+"\t"+ stripped[6]+"\t"+TE)
        elif line.startswith("Kimura"):
                stripped=line.split("=")
                if float(stripped[1]) <= float(arg.k):
                    kimuraval = stripped[1].rstrip("\n")
                    name[Id]=kimuraval
        elif line.startswith("Transitions"):
                strippedt=line.split(" ")
                trans=strippedt[4]
#                print (trans)
                Tranver[Id]=trans
                
oWriter=open(arg.o, "w")
oWriter.write("Scaffold"+"\t"+"beg"+"\t"+"end"+"\t"+"type"+"\t"+"kimuraval"+"Transitions/transversions"+"\n")
bWriter=open(arg.b, "w")
for rec in name:
        #print (rec+"\t" +name[rec]+"\t" + Tranver[rec])
        oWriter.write(rec+"\t" +name[rec]+"\t" + Tranver[rec]+"\n")
        bWriter.write(rec+"\t" +name[rec]+"\t" + Tranver[rec]+"\n")
        #print(rec +name[rec])
