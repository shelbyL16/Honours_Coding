from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime
startTime = datetime.now()
print("It started at this time",startTime)

input_file="C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\sequence2.fasta"
output_file = "C:\\Users\\slabu\\Desktop\\University\\Work\\Project\\14 August Starting fresh\\res.fasta"
fasta_sequences = SeqIO.parse(open(input_file),'fasta')

i = 0
FRAME_SIZE = 7
sequence_description_map = dict()
MAX=6000 #this is the number of sequences you want to read in from the file if you don't want to read in the whole file
#DO NOT MAKE 1 AS YOU WILL BE STUCK IN AN INFINITE LOOP
ratio=0.8#this is the percentage of sequences you will keep at the end of each run
target=400#this is the number you would like the sequences to be less than


#This is the method that loads the data from the file to the sequence map for manipulation
def load_data():
    sequence_map = dict()#this is the dictionary with the sequence IDs, we have two different dictionaries!
    sequence_size_map = dict()#computing the length, then if have length, store the id to the length, keeps similar sequence lengths together
    #not using this in any way I think
    with open(output_file) as out_file:
        i = 0
        for fasta in fasta_sequences:#this will loop through the file
            if i == MAX:#this is used if you don't want to read in the whole file
                break
            name, sequence, description = fasta.id, str(fasta.seq), fasta.description#getting the details from the file
            #the name is the sequenceID
            #sequence is the dna sequence
            #description is the strain and where it was found I think
            sequence_map [name] = sequence#map the sequences to their names in the dictionary
            sequence_description_map[name]=description#save the desc to the ID, used later when saving to the file
            #this can be used if want to delete sequences first according to their lenght
            length = len(sequence)#stored length as a variable so that we can take away according to the lengths at some point
            if sequence_size_map.get(length):#fetch it if it is there, incase it aint there guys, do the false condition
                sequence_size_map[length].append(name)#adding the name to the list of sequences with the same length
            else: 
                sequence_size_map[length] = [name]#this is the first one in this case so add it to the dictionary with the length as the key
            sequence_map [name] = sequence#map the sequences to their names
            i+=1
        #for key in sequence_size_map:
        # print(key)
        return sequence_map

#this is the method that calculates the score for the sequence similarity
def get_sequence_similarity(sequence1, sequence2):
    #saving subsequences found in a dictionary so that if we encounter it again at a later stage we can skip it
    #increase efficiency of the program as we would have already done the work
    #this is for the possible repitition of the sequences that we have accounted forS
    score = 0
    count = 0
    processed_sub_sequence = dict()#this is the dictionary that is used to keep track of the sequences repitition
    for i in range(len(sequence1) - FRAME_SIZE):#default is zero then it will go to the end of the sequences, looping through the first sequence
        nt_sequence = sequence1[i:i+FRAME_SIZE]#substring from i to the frame size to keep the matches to 7 bases as suggested
        if processed_sub_sequence.get(nt_sequence):#If we have previously found this sequence, don't do it again
            count = processed_sub_sequence[nt_sequence]#count equals the number at this one for the sequence
        else:#havent found it previously
            count = sequence2.count(nt_sequence)#count how many times it is there
            processed_sub_sequence[nt_sequence] = count#store it to the map
            
        score += count#increase the score with count
    #return score
    return (score/(len(sequence1)+len(sequence2)))#dividing by the lengths to "normalise" in relation to big and small sequences
#complexity = 100
#seq1 = "AGGTAGGATTCATT"
#seq2 = "AGGTAGGATCCATT"
#seq1 = sequence_map[keys[0]]#[:10000]
#seq2 = sequence_map[keys[1]]#[:10000]
#print(seq1)
#print(seq2)

#print(get_sequence_similarity(seq1, seq2))
#print(len(sequence_map[keys[0]]))
#print(len(sequence_map[keys[1]]))

#comparison uses get_sequence_similarity to calculate similarity between the sequences
def comparison(keys,sequence_map):

    score_matrix = dict()#score matrix to store the score answers
    seq1 = sequence_map[keys[0]]#get the first sequence from the sequence mapped parsed as we always compare to the first one
    for i in range(1,len(keys)):#go from the next sequence to the length of keys, keys is a shortened list
        seq2 = sequence_map[keys[i]]#get the relevant sequence to compare sequence 1 to in the for loop
        score_matrix [(keys[0],keys[i])] = get_sequence_similarity(seq1, seq2)#store the value of the score for that sequence at this dictionary entry
    
    sorted_score = sorted(score_matrix .keys(), key=lambda x: score_matrix [x])#sort the sequences to get the most dissimilar ones at the top
    sorted_keys = [x[1] for x in sorted_score]#adding the first sequence back at the second position so we can compre it - increase accuracy
    return [sorted_keys[0], keys[0]] + sorted_keys[1:]#improving accuracy
    # return sorted_score

# def list_scores(sc):#take in the score matrix for this function getting the average of the scores
#     score_similarity = dict()
#     temp = 0
#     for i in range(MAX):#go through sequence map to get the sequence IDs
#          for j in sc.keys():#got to loop through the matrix to find the corresponding sequence ID scores
#             if keys[i]  in j:
#                 temp = temp+sc[j]
#          score_similarity[keys[i]] = temp
#          temp = 0
#     sorted_score = sorted(score_similarity.keys(), key=lambda x: score_similarity[x])
#     return sorted_score

#this writes the sequences we have found to the file at the end
def sequences_write_to_file(num,sorted_keys,sequence_map):
    temp = ""
    records = []
    for i in range(len(sorted_keys)):#loop through the sorted keys
        temp = sequence_map[sorted_keys[i]]#find the entry from the sorted keys that maps to the sequence we need from sequence_map
       # seq_record.seq = Seq(re.sub('[^GATC]',"",str(sequence).upper()))
        record = SeqRecord(Seq(temp),  str(sorted_keys[i]), description=sequence_description_map[sorted_keys[i]])#create a record, this is needed to write to a fasta file
        records.append(record)#append the record to the array
    SeqIO.write(records, output_file, "fasta")#write the records to the fasta file

#need to update our sequence map as we go when we sort the keys
def update_sequence_map(sorted_keys,length,sequence_map):#length is sequences you want and sequence map
    res = dict()
    for i in range(length):#loop through according to length
        res[sorted_keys[i]]=sequence_map[sorted_keys[i]]#add to res dictionary what was at sequence map
    return res

#call this method to make the program run
def main(length):#this is the target length of sequences
    sequence_map = load_data()#call load data first as need to load it into the sequence map
    keys = list(sequence_map.keys())#need to do it first time to see it, get the keys as a list so we can manipulate them
    sorted_keys = comparison(keys, sequence_map)#comparison uses get_sequence_similarity to calculate similarity between the sequences and return a sorted list of keys
    sequence_map = update_sequence_map(sorted_keys, int(len(sequence_map)*ratio),sequence_map)#call update data to reload the data losig the ones that are similar according to your ratio inputted
    print("this is the length", len(sequence_map), "time", (datetime.now() - startTime))#statements just to check progress and timing


    while len(sequence_map)>length:#person inputs the number of sequences they want and the map is still longer than wanted
        keys = list(sequence_map.keys())#list the keys again
        sorted_keys = comparison(keys, sequence_map)#sort the keys, comparison returns sorted keys with get_similarity method
        sequence_map = update_sequence_map(sorted_keys, int(len(sequence_map)*ratio),sequence_map)#update sequence data according to the ration you entered
        print("this is the length", len(sequence_map), "time", (datetime.now() - startTime))#statements just to check progress and timing
    keys = list(sequence_map.keys())#list the keys again form the new sequence map we have calculated from above
    sorted_keys = comparison(keys, sequence_map)#sort the keys, comparison returns sorted keys with get_similarity method
    sequences_write_to_file(length,sorted_keys[:len(sequence_map)],sequence_map)#write to file when we have the sequences we wanted

main(target)#this is the method call that starts the whole program
print("It took:", (datetime.now() - startTime))