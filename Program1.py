from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

input_file="C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\sequence.fasta"
output_file = "C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\res.fasta"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')

i = 0
FRAME_SIZE = 7
sequence_description_map = dict()

#This is the method that loads the data from the file to the sequence map
def load_data():
    sequence_map = dict()#this is the dictionary with the sequence IDs, we have two different dictionaries!
    sequence_size_map = dict()#computing the length, then if have length, store the id to the length, keeps similar sequence lengths together
    with open(output_file) as out_file:
        for fasta in fasta_sequences:
            name, sequence, description = fasta.id, str(fasta.seq), fasta.description
            #the name is the sequenceID
            #sequence is the dna sequence
            sequence_map [name] = sequence#map the sequences to their names in the dictionary
            sequence_description_map[name]=description#save the desc to the ID

            length = len(sequence)#stored length as a variable so that we can take away according to the lengths at some point
            if sequence_size_map.get(length):#fetch it if it is there, incase it aint there guys, do the false condition
                sequence_size_map[length].append(name)#adding the name to the list of sequences with the same length
            else: 
                sequence_size_map[length] = [name]#this is the first one in this case so add it to the dictionary with the length as the key
            sequence_map [name] = sequence#map the sequences to their names
        
        #for key in sequence_size_map:
        # print(key)
        return sequence_map

def get_sequence_similarity(sequence1, sequence2):
    #saving subsequences found in a dictionary so that if we encounter it again at a later stage we can skip it
    #increase efficiency of the program as we would have already done the work
    #this is for the possible repitition of the sequences that we have accounted forS
    score = 0
    processed_sub_sequence = dict()#this is the dictionary that is used to keep track of the sequences repitition
    for i in range(len(sequence1) - FRAME_SIZE):#default is zero then it will go to the end of the sequences, looping through the first sequence
        nt_sequence = sequence1[i:i+FRAME_SIZE]#substring from i to the frame size to keep the matches to 7 bases as suggested
        j = 0
        # if processed_sub_sequence.get(nt_sequence)  == True:#If we have previously found this sequence, don't do it again
        #     continue
        processed_sub_sequence[nt_sequence] = True#if it is not found then the score is zero
        #this is the one where it will continue to go with the nucleotides if it has found a match and sees how long they are
        while j != -1: #will terminate when it is not found
            j = sequence2.find(nt_sequence, j) #find the sequence in the second sequence
            if j != -1 :
                score += 1#increase the score cause we have found the matches
                j += 1
    score = (score/(len(sequence1)+len(sequence2)))*100
    return score

sequence_map = load_data()
keys = list(sequence_map.keys())
# seq1 = "AGGTAGGATTCATT"
# seq2 = "AGGTAGGATCCATT"
#seq1 = sequence_map[keys[0]]#[:10000]
#seq2 = sequence_map[keys[1]]#[:10000]
#print(seq1)
#print(seq2)

#print(get_sequence_similarity(seq1, seq2))
#print(len(sequence_map[keys[0]]))
#print(len(sequence_map[keys[1]]))

def comparison():

    score_matrix = dict()#score matrix
    for i in range(5):
        for j in range(i+1,5):#two inner loops to go through the sequences and compare them
            seq1 = sequence_map[keys[i]]#get the relevant sequences
            seq2 = sequence_map[keys[j]]#get the relevant sequences
            score_matrix [(keys[i],keys[j])] = get_sequence_similarity(seq1, seq2)#store the value of the score for that sequence at this dictionary
    return score_matrix#return the score matrix for all the compared sequences

def list_scores(sc):#take in the score matrix for this function getting the average of the scores
    score_similarity = dict()
    temp = 0
    for i in range(5):#go through sequence map to get the sequence IDs
         for j in sc.keys():#got to loop through the matrix to find the corresponding sequence ID scores
            if keys[i]  in j:
                temp = temp+sc[j]
         score_similarity[keys[i]] = temp
         temp = 0
    sorted_score = sorted(score_similarity.keys(), key=lambda x: score_similarity[x])
    return sorted_score
   
def sequences_write_to_file(num,sorted_keys):
    temp = ""
    records = []
    for i in range(num):#loop through to find the sequences that match to that key
        temp = sequence_map[sorted_keys[i]]
       # seq_record.seq = Seq(re.sub('[^GATC]',"",str(sequence).upper()))
        record = SeqRecord(Seq(temp),  str(sorted_keys[i]), description=sequence_description_map[sorted_keys[i]])
        records.append(record)
    SeqIO.write(records, output_file, "fasta")#write it to the fasta file

# comparison()
sequences_write_to_file(3,list_scores(comparison()))