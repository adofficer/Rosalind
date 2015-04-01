import math
import timeit

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

generations = [1,1]
def mortal_fib(i, j):
    count = 2
    while (count < i):
        if (count < j):
            generations.append(generations[-2] + generations[-1]) #recurrence relation before rabbits start dying (simply fib seq Fn = Fn-2 + Fn-1)
        elif (count == j or count == j+1):
            generations.append((generations[-2] + generations[-1]) - 1)#Fn = Fn-2 + Fn-1 - 1
        else:
            generations.append((generations[-2] + generations[-1]) - (generations[-(j+1)])) #Our recurrence relation here is Fn-2 + Fn-1 - Fn-(j+1)
        count += 1
    return (generations[-1])

def expect_offspring(AA_AA, AA_Aa, AA_aa, Aa_Aa, Aa_aa, aa_aa):
    total = (AA_AA+AA_Aa+AA_aa+Aa_Aa+Aa_aa+aa_aa)
    return 2*(1*AA_AA + 1*AA_Aa + 0.75*Aa_Aa + 1*AA_aa + 0.5*Aa_aa)


def mrna(prot_string):
    possib_codons = { "I":3, "Y":2, "F":2, "L":6, "M":1, "V":4, "S":6, "P":4, "T":4, "A":4, "H":2, "Q":2, "N":2, "K":2, "D":2, "E":2, "C":2, "W":1, "R":6, "G":4 }
    possib_num = 1
    for char in str(prot_string):
        possib_num = (possib_num*possib_codons[char] % 1000000)
    print (possib_num * 3 % 1000000)

def enum_gene_order(n):
    if n>1:
        permute_list = enum_gene_order(n-1)
        next_permute_list = []
        for prior_perm in permute_list:
            for i in range(0, len(prior_perm)+1):
                new_entry = list(prior_perm)
                new_entry.insert(i, n)
                next_permute_list.append(new_entry)
        permute_list = list(next_permute_list)
        return permute_list
    else:
        return [[1]]

def print_nicely(input_list):
    for entry in input_list:
        for i in range(0, len(entry)):
            print entry[i],
        print

def cons(dna_dict):
    for key in dna_dict:
        len_list = len(dna_dict[key])
        break
    a_list = [0 for x in range(0, len_list)]
    c_list = [0 for x in range(0, len_list)]
    g_list = [0 for x in range(0, len_list)]
    t_list = [0 for x in range(0, len_list)]
    for i in range(0, len_list):
        for key in dna_dict:
            if dna_dict[key][i] == "T":
                t_list[i] += 1
            if dna_dict[key][i] == "C":
                c_list[i] += 1
            if dna_dict[key][i] == "G":
                g_list[i] += 1
            if dna_dict[key][i] == "A":
                a_list[i] += 1
    print a_list
    print c_list
    print g_list
    print t_list


cons(fasta_parse(....))
