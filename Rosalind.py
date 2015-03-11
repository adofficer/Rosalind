def gc_cont():
    file = open(raw_input("File name? "), "r")
    dummy = file.read()
    file.close()
    input_list = dummy.splitlines()
    gc_count = 0.00000
    str_len = 1.00000
    name = ""
    best_name = ""
    best_gc = 0.00000
    for item in input_list:
        if item[0] == ">":
            if (gc_count/str_len) > best_gc:
                best_name = name
                best_gc = (gc_count/str_len)
            name = item[0:]
            gc_count = 0.00000
            str_len = 0.00000
        else:
            gc_count += (item.count("C") + item.count("G"))
            str_len += (item.count("C") + item.count("G") + item.count("A") + item.count("T"))
    if (gc_count/str_len) > best_gc:
        best_name = name
        best_gc = (gc_count/str_len)
    print best_name[1:] + " " + str(100*best_gc)[0:9]
    
def point_mut():
    import itertools
    file = open(raw_input("File name? "), "r")
    dummy = file.read()
    input_list = dummy.splitlines()
    ham_dist = 0
    for item1, item2 in zip(input_list[0], input_list[1]):
        if item1 != item2:
            ham_dist += 1
    print ham_dist
    file.close()


def mendel1():
    file = open(raw_input("File name? "), "r")
    dummy = file.read()
    kmn_list = dummy.split(" ")
    kmn_tot = float(kmn_list[0]) + float(kmn_list[1]) + float(kmn_list[2])
    for i in range(0, len(kmn_list)):
        kmn_list[i] = float(kmn_list[i])
    key_list = ["het v het", "het v homorec", "homorec v homorec"]
    sum = 0
    prob_rec_list = [0.25, 0.5, 1]
    prob_mate_list = [kmn_list[1]/kmn_tot*(kmn_list[1]-1)/(kmn_tot-1),
                      2*kmn_list[1]/kmn_tot*kmn_list[2]/(kmn_tot-1),
                      kmn_list[2]/kmn_tot*(kmn_list[2]-1)/(kmn_tot-1)]
    for i in range(0, len(prob_rec_list)):
        sum += prob_rec_list[i] * prob_mate_list[i]
    print 1-sum
    file.close()

rna_table = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def rnatrans():
    file = open(raw_input("File name? "), "r")
    rna_seq = file.read()
    file.close()
    prot_seq = ""
    for i in range(0, len(rna_seq)/3):
        if rna_table[rna_seq[3*i:3*i+3]] == "STOP":
            break
        else:
            prot_seq = prot_seq + rna_table[rna_seq[3*i:3*i+3]]
    print prot_seq.upper()

def dna_subs():
    file = open(raw_input("File name? "), "r")
    dummy = file.read().splitlines()
    file.close()
    dna_ref_seq = dummy[0]
    probe_seq = dummy[1]
    location_matches = []
    for i in range(1, len(dna_ref_seq)-3):
        if probe_seq[0:4] == dna_ref_seq[int(i):int(i)+4]:
            location_matches.append(i)
    print location_matches
    for item in location_matches:
        if probe_seq == dna_ref_seq[item:item+len(probe_seq)]:
            print item+1,

gc_cont()
