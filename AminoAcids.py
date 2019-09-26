from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime
from  collections import Counter
import csv



ans = ""
input_file="C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\FUBAR\\2ndSeptember\\5thSeptember\\13th September\\ShelbyAnalysis\\NewData\\AminoAcidsHERE.fas"
out_file = "C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\FUBAR\\2ndSeptember\\5thSeptember\\13th September\\ShelbyAnalysis\\NewData\\AminoAcidProfiles.txt"
csv_file = "C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\FUBAR\\2ndSeptember\\5thSeptember\\13th September\\ShelbyAnalysis\\NewData\\AminoAcidProfilesCSV.csv"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
map_name_seq = dict()
map_desc_seq = dict()
places_dict = dict()
amino_acids=[]
ans = []

i=0
with open(input_file) as input_file:
    for fasta in fasta_sequences:
        #print("this is the id",fasta.id[0:10])
        name, sequence,desc = fasta.id[0:10], str(fasta.seq),fasta.id[10::]
        map_name_seq [name] = sequence#map the sequences to their names in the dictionary
        map_desc_seq [name] = desc#map the desc to their names in the dictionary
        amino_acids.append(sequence)
        
temp = ""
for i in range (0,679):#get an amino acid sequence- amino acids are 679 in len, but there are 99
    for j in range (0,99):#this is the number of amino acids
        seq = amino_acids[j]#this is the amino acid sequence
        temp = temp+seq[i]
    ans.append(temp)
    temp = ""


#print(len(ans))
#print(len(ans[0]))

d = {}
ans_arr={}
count = 1

for temp in ans:#count occurences of the string
    #print(count,temp)#know it has all the amino acids at that site
    for t in temp:
        d[t] = temp.count(t)
    ans_arr[count]=d
    count+=1
    d={}


AminoAcidDictionary = {
'A' : "Alanine",
'R' : "Arginine",
'N' : "Asparagine",
'D' : "Aspartic Acid",
'B' : "Asparagine Aspartic Acid",
'C' : "Cysteine",
'E' : "Glutamic Acid",
'Q' : "Glutamine", 
'Z' : "Glutamine Glutamic acid", 
'G' : "Glycine",
'H' : "Histidine", 
'I' : "Isoleucine",
'L' : "Leucine",
'K' : "Lysine",
'M' : "Methionine",
'F' : "Phenylalanine",
'P' : "Proline",
'S' : "Serine",
'T' : "Threonine",
'W' : "Tryptophan", 
'Y' : "Tyrosine", 
'V' : "Valine",
'-' : "Space"}

f= open(out_file,"w+")
y = ""
c= ""
#want to write a return method that will write what we have to a file and loop up the amino acids  
for i in range(1,len(ans_arr)):#loop through the big array
    strt = ""
    temp_d = {}#need a temp dict
    temp_d = ans_arr[i]#save the dict at i to this temp one
    for t in temp_d:#need to loop through the dictionary 
    #t is the amino acid so thats cool but now I need to get the number associated with it 
        #print(temp_d[t])
        strt = strt +  AminoAcidDictionary.get(t) + " : " + str(temp_d[t])  +"       "   #trying to get it to name the amino acid
    y=y+ "Site" + str(i) +"                 "+ strt+ "\n"

f.write(str(y))       
f.close()

#####Trying to write to a fasta file 
'''
csvData = [['Person', 'Age'], ['Peter', '22'], ['Jasmine', '21'], ['Sam', '24']]
with open(csv_file, 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(csvData)
csvFile.close()
'''