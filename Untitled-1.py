from Bio import SeqIO
input_file="C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\sequence.fasta"
output_file = "C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\res.fasta"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')

i = 0
def load_data():
    sequence_map = dict()
    sequence_size_map = dict()#computing the length, then if have length, store the id to the length
    with open(output_file) as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            #the name is the sequenceID
            #sequence is the dna sequence
            # new_sequence = some_function(sequence)
            # write_fasta(out_file)
            sequence_map [name] = sequence#map the sequences to their names
            length = len(sequence)#stored length as a variable
            if sequence_size_map.get(length):#fetch it if it is there, incase it aint there guys
                sequence_size_map[length].append(name)#adding the name to the list of sequences with the same length
            else: 
                sequence_size_map[length] = [name]#this is the first one in this case
            sequence_map [name] = sequence#map the sequences to their names
            #print(name, i)
            #print(len(sequence),sequence, "\n")
            #i += 1
        
        #for key in sequence_size_map:
        # print(key)