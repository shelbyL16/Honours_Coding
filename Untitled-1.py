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
            length = len(sequence)#stored length as a variable so that we can take away according to the lengths at some point
            if sequence_size_map.get(length):#fetch it if it is there, incase it aint there guys
                sequence_size_map[length].append(name)#adding the name to the list of sequences with the same length
            else: 
                sequence_size_map[length] = [name]#this is the first one in this case so add it to the dictionary
            sequence_map [name] = sequence#map the sequences to their names
            #print(name, i)
            #print(len(sequence),sequence, "\n")
            #i += 1
        
        #for key in sequence_size_map:
        # print(key)
        return sequence_map
FRAME_SIZE = 3

def get_sequence_similarity(sequence1, sequence2):
    #saving subsequences found in a dictionary so that if we encounter it again at a later stage we can skip it
    #increase efficiency of the program as we would have already done the work
    #this is for the possible repitition of the sequences
    score = 0
    processed_sub_sequence = dict()
    for i in range(len(sequence1) - FRAME_SIZE):#default is zero then it will go to the end of the sequenceS
        nt_sequence = sequence1[i:i+FRAME_SIZE]#substring from i to the frame size
        j = 0
        if processed_sub_sequence.get(nt_sequence)  == True:
            continue
        processed_sub_sequence[nt_sequence] = True#if it is not found then the score is zero
        while j != -1: #will terminate when it is not found
            j = sequence2.find(nt_sequence, j) #find the sequence in the second sequence
            if j != -1 :
                score += 1
                j += 1
    return score

sequence_map = load_data()
keys = list(sequence_map.keys())
# seq1 = "AGGTAGGATTCATT"
# seq2 = "AGGTAGGATCCATT"
seq1 = keys[0]
seq2 = keys[1]
print(sequence_map[seq1])
print(sequence_map[seq2])

print(get_sequence_similarity(seq1, seq2))