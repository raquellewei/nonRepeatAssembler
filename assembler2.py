#test if this is a "good-enough" version
#1. filter with filter functions 
#2. filter repeats
#3. concatenate
#4. filter again with global filter
#5. filter again to get rid of additional repeats

from nrpFinder.finder import nrp_finder as nrpfinder
from nrpFinder import utils
from nrpFinder import hgraph, recovery
import re
from random import choice, randint, shuffle
import numpy as np
import csv
import networkx as nx 

#parse csv files and return dict
def parse_csv(file):
	parts_dict = {}
	with open(file, "r") as infile:
		reader = csv.reader(infile)
		for row in reader:
			if not row[0].isalpha():
				parts_dict[row[0]] = row[1]
	return parts_dict

def csv_to_fa(csv_file):
	fa_file_name = csv_file.replace('csv', 'fa')
	with open(fa_file_name, "w") as fa_file:
		with open(csv_file, "r") as infile:
			infile.readline()
			for line in infile:
				index,seq = line.strip().split(',')
				fa_file.write('> '+index+'\n')
				fa_file.write(seq+'\n')
	return fa_file_name


#some filter functions
def gc_content(parts_list, percent, sign):
	parts_copy = parts_list[:]
	for seq in parts_copy:
		gc_count = 0
		for n in seq:
			if n.lower() == "g" or n.lower() == "c":
				gc_count+=1
		gc_content = float(gc_count)/len(seq)
		if sign == ">":
			if gc_content < percent:
				parts_copy.remove(seq)
		elif sign == "<":
			if gc_content > percent:
				parts_copy.remove(seq)
		else:
			print "sign has to be either '<' or '>'"
	return parts_copy


def does_not_allow(parts_list, *motifs):
	parts_copy = parts_list[:]
	for seq in parts_copy:
		for motif in motifs:
			if re.match(motif, seq, re.IGNORECASE):
				parts_copy.remove(seq)
	return parts_copy

def check_motif(str, motif):
	#returns if this motif is found 
	if re.match(motif, str, re.IGNORECASE):
		return True
	else:
		return False


#get non-repetitive parts

#assemble
def assemble_naive(order, homology, *motifs):
	id_dict, seq_list = setup(*motifs)
	homology_graph = hgraph.get_homology_graph(seq_list[:], [], homology, False, 'seq_file.txt', verbose = False)
	indi_seqs_indx = recovery.get_recovered_non_homologs(homology_graph, 'graph_file.txt', None, verbose = False)
	updated_dict = standardize_lists(id_dict, indi_seqs_indx)
	#assemble
	order = order.split(',')
	assembled = []
	while True:
		seq = ""
		for number in order:
			if len(updated_dict[number]) == 0: #toolbox empty
				filtered = does_not_allow(assembled, *motifs)
				return filtered
			else:
				seq_id = updated_dict[number].pop()
				seq += seq_list[seq_id]
		assembled.append(seq)
	print "raw assemble: {}".format(len(assembled))
	filtered = does_not_allow(assembled, *motifs)
	return filtered


def assemble_max_utilize(order, homology, *motifs):
	#build homology graph, check motif at every step, maximize usage of each part
	id_dict, seq_list = setup(*motifs)
	G = hgraph.get_homology_graph(seq_list[:], [], homology, False, 'seq_file.txt', True)
	updated_dict = standardize_lists(id_dict, G.nodes())
	
	G = assign_attribute(G,updated_dict)

	order = order.split(',')
	assembled = []
	while True:
		if len(updated_dict[order[0]]) == 0:
			return assembled
		else:
			i = choice(updated_dict[order[0]])
			seq_prev = seq_list[i]
			for j in list(G.neighbors(i)):
				where = G.nodes[j]['type']
				G.remove_node(j)
				updated_dict[where].remove(j)
			G.remove_node(i)
			updated_dict[order[0]].remove(i)
		for t in order[1:]:
			if len(updated_dict[t]) == 0:
				return assembled
			else:
				used = []
				while not len(updated_dict[t]) == len(used):
					seq_id = choice(updated_dict[t])
					while seq_id in used:
						seq_id = choice(updated_dict[t])
					#concatenate and test
					seq_test = seq_prev + seq_list[seq_id]
					for motif in motifs:
						if check_motif(seq_test, motif): #if motif is found
							print "found motif {}".format(motif)
							print len(used)
							used.append(seq_id)
							break
					else:
						#concatenate. remove the rest
						seq_prev = seq_test
						for n in list(G.neighbors(seq_id)):
							where = G.nodes[n]['type']
							G.remove_node(n)
							updated_dict[where].remove(n)
						G.remove_node(seq_id)
						updated_dict[t].remove(seq_id)
						used = []
						break
				else:
					return assembled
		assembled.append(seq_prev)
	return assembled
					


def assign_attribute(h_graph, id_dict):
	for key in id_dict:
		for i in id_dict[key]:
			if i in h_graph.nodes():
				h_graph.nodes[i]['type'] = key
	return h_graph



def standardize_lists(id_dict, no_repeat_list):
	for key in id_dict:
		id_dict[key] = filter(lambda i : i in no_repeat_list, id_dict[key])
	return id_dict
#main
def main():
	pass


if __name__ == '__main__':
    main() 


#test functions

def setup(*motifs):
	type1_fa = csv_to_fa("NonRepetitiveParts_Type1.csv")
	type2_fa = csv_to_fa("NonRepetitiveParts_Type2.csv")
	type3_fa = csv_to_fa("NonRepetitiveParts_Type3.csv")

	seq_list1 = utils.get_fasta_seq_list(type1_fa)
	seq_list2 = utils.get_fasta_seq_list(type2_fa)
	seq_list3 = utils.get_fasta_seq_list(type3_fa)
	seq_list1 = does_not_allow(seq_list1, *motifs)
	seq_list2 = does_not_allow(seq_list2, *motifs)
	seq_list3 = does_not_allow(seq_list3, *motifs)
	seq_list = seq_list1+seq_list2+seq_list3

	seq_1_id = range(len(seq_list1))
	seq_2_id = range(len(seq_list1), len(seq_list1)+len(seq_list2))
	seq_3_id = range(len(seq_list1)+len(seq_list2), len(seq_list1)+len(seq_list2)+len(seq_list3))
	id_dict = {'1':seq_1_id, '2': seq_2_id, '3':seq_3_id}
	return id_dict, seq_list


####
max_assembled = assemble_max_utilize('2,3,2,1', 15, 'ATT')
assembled = assemble_naive('2,3,2,1', 15, 'ATT')

print "max assemble yields: {}".format(len(max_assembled))
print "naive yields: {}".format(len(assembled))

non_homologs_max = nrpfinder(max_assembled,[], 15, verbose=False)
non_homologs = nrpfinder(assembled,[], 15, verbose=False)

print len(non_homologs_max)
print len(non_homologs)

#print non_homologs_max
#print non_homologs
#indi_seqs_indx = recovery.get_recovered_non_homologs(homology_graph, 'graph_file.txt', None, verbose = False)
#id_dict_updated = standardize_lists(id_dict, indi_seqs_indx)
#assemble_using_graph('1,2,3', 13)
#print assembled
#print non_homologs



#indi_seqs_indx = recovery.get_recovered_non_homologs(homology_graph, 'graph_file.txt', vercov_func=None, verbose=True)
#print indi_seqs_indx
#print seq_list
#filtered_list = does_not_allow(seq_list, '^Gggc')
#print len(seq_list)
#print len(filtered_list)
