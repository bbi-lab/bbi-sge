#!/usr/bin/env python

# cigar

# amplicon

# exp

def import_dict(my_dict_file):
	my_dict = {}
	dict_file = open(my_dict_file, 'r')
	for line in dict_file:
		data = line.strip().split('\t')
		my_dict[data[0]] = int(data[1])
	return my_dict

#make a dictionary from a matrix file in which each key is defined from the first tab, and it is paired with a list consisting of tabs 2 through N
#will also spit out a header line first in a two part tuple if there is a header in the dictionary (will assume there is)
def import_matrix2dict(my_matrix_file, header=True):
	my_dict = {}
	my_header = []
	my_matrix = open(my_matrix_file, 'r')
	#does the first line get put into the dictionary or the header?
	if header:
		line_count = 0
	else:
		line_count = 1

	for line in my_matrix:
		if line_count == 0:
			#header returns as list
			header_list = line.strip().split('\t')
			line_count +=1
		else:
			data = line.strip().split('\t')	#dictionary now takes form variant:list of counts/features
			my_dict[data[0]] = data[1:]
	if header:
		return (header_list,my_dict)
	else:
		return my_dict


#functions to merge dictionaries, normalize dictionaries, take ratios between counts in each dictionary, and add CADD_scores
#give two dictionaries, and choose if you want to merge on "A"'s data, the "intersect", or the "union":

#for now assumes that each dict has counts + some other information about the variant (cigar string, annotations, etc.)
def merge_dict_counts(dict_a, dict_b):
	merged_dict = {}
	dict_a_keys = dict_a.keys()
	for a_key in dict_a_keys:
		if a_key in dict_b:
			merged_dict[a_key] = [dict_a[a_key][0],dict_b[a_key][0]]
		else:
			merged_dict[a_key] = [dict_a[a_key][0],0]
	for b_key in dict_b:
		if b_key not in merged_dict:
			merged_dict[b_key] = [0,dict_b[b_key][0]]
	return merged_dict

def sum_dicts(dict_a, dict_b, merge_type = 'A'):
	merged_dict = {}
	if merge_type == 'A':
		dict_a_keys = dict_a.keys()
		for a_key in dict_a_keys:
			if a_key in dict_b:
				merged_dict[a_key] = dict_a[a_key]+dict_b[a_key]
			else:
				merged_dict[a_key] = dict_a[a_key]
		return merged_dict
	elif merge_type == 'union':
		dict_a_keys = dict_a.keys()
		dict_b_keys = dict_b.keys()
		for a_key in dict_a_keys:
			if a_key in dict_b:
				merged_dict[a_key] = dict_a[a_key]+dict_b[a_key]
			else:
				merged_dict[a_key] = dict_a[a_key]
		for b_key in dict_b_keys:
			if b_key in merged_dict:
				continue
			else:
				merged_dict[b_key] = dict_b[b_key]
		return merged_dict
	elif merge_type == 'intersect':
		dict_a_keys = dict_a.keys()
		for a_key in dict_a_keys:
			if a_key in dict_b:
				merged_dict[a_key] = dict_a[a_key]+dict_b[a_key]
		return merged_dict
	else:
		print "Please specify a valid merge type, either 'A', 'union', or 'intersect'."

#to divide dict_b values by dict_a values -- 'A' will perform for all values for A condition, "union" will pseudocount +1 any missing values and perform,
# "intersect" will perform only for all values represented in both -- threshold is the minimum number of reads required in dict_a to include in output
def ratio_dicts(dict_a, dict_b, merge_type = 'A', threshold = 6):
	merged_dict = {}
	if merge_type == 'A':
		dict_a_keys = dict_a.keys()
		for a_key in dict_a_keys:
			if dict_a[a_key] >= threshold:
				if a_key in dict_b:
					merged_dict[a_key] = float(dict_b[a_key])/dict_a[a_key] 
				else:
					merged_dict[a_key] = 0
		return merged_dict
	elif merge_type == 'union':
		dict_a_keys = dict_a.keys()
		dict_b_keys = dict_b.keys()
		for a_key in dict_a_keys:
			if dict_a[a_key] >= threshold:
				if a_key in dict_b:
					merged_dict[a_key] = dict_b[a_key]/dict_a[a_key]
				else:
					merged_dict[a_key] = 0
		for b_key in dict_b_keys:
			if dict_b[b_key] >= threshold:
				if b_key in merged_dict:
					pass
				else:
					#if not in dict_a, just report value in dict_b/1
					merged_dict[b_key] = dict_b[b_key]/1
		return merged_dict
	elif merge_type == 'intersect':
		dict_a_keys = dict_a.keys()
		for a_key in dict_a_keys:
			if a_key in dict_b and dict_a[a_key] >= threshold:
				merged_dict[a_key] = dict_b[a_key]/dict_a[a_key]
		return merged_dict
	else:
		print "Please specify a valid merge type, either 'A', 'union', or 'intersect'."
		
#normalize_dict -- will divide each value by the sum of all values to normalize for sequencing coverage
def normalize_dict(dict_a):
	normalized_dict = {}
	dict_a_reads = 0
	dict_a_keys = dict_a.keys()
	dict_a_keys.sort()
	for a_key in dict_a_keys:
		dict_a_reads += dict_a[a_key]
	for a_key in dict_a_keys:
		normalized_dict[a_key] = float(dict_a[a_key])/dict_a_reads
	return normalized_dict

#needs full path for cadd_score_file, chr_coord should be the start coordinate of the wt sequence ##change this post editing annotate script
def annotate_w_cadd(my_dict, cadd_score_file, wt_seq, chr_coord, dict_header=True, cadd_header=True):

	index_0 = int(chr_coord)
	#format cadd_data into a lookup table
	open_cadd_scores = open(cadd_score_file, 'r')
	cadd_data = {}

	line_count = 0
	if cadd_header != True:
		line_count = 2
	for line in open_cadd_scores:
		line_count += 1
		if line_count > 2:
			variant_data = my_line.strip().split('\t')
			variant_key = str(variant_data[1])+str(variant_data[4])
			#maps each variant to all of its cadd data in a list format.
			if variant_key not in cadd_data:
				cadd_data[variant_key] = variant_data
			
	#go through read_dict and find reads of the wt length
	annotated_dict = {}
	full_length_seqs = 0
	single_nuc_subs = 0
	wt_len = len(wt_seq)
	SNV_in = 0
	SNV_out = 0
	snv_keys_not_in_CADD = []
	for seq in my_dict:
		if len(seq)==len(wt_seq):
			#count substitutions from wt seq
			full_length_seqs += 1
			subs = 0
			sub_index = 0
			for i in range(0,wt_len):  		#aaron recommends itertools' zip module
				if wt_seq[i] != seq[i]:
					subs+=1
					sub_index = i
					if subs>1:
						sub_index = 0
						continue
			if subs>1:
				continue
			elif subs==0:
				annotated_dict[seq] = my_dict[seq] + ['wt','wt','wt','wt','wt','wt','wt','wt','wt','wt']
			elif subs==1:
				#print seq
				single_nuc_subs += 1
				position = index_0+sub_index 
				variant = seq[sub_index]
				snv_key = str(position)+str(variant)
				if snv_key in cadd_data:
					SNV_in += 1
					variant_CADD_data = cadd_data[snv_key]
					selected_variant_CADD_data = [variant_CADD_data[0],variant_CADD_data[1],variant_CADD_data[2],variant_CADD_data[4],variant_CADD_data[10],variant_CADD_data[12],variant_CADD_data[110],variant_CADD_data[111],variant_CADD_data[114],variant_CADD_data[115]]
					prev_data_list = my_dict[seq]
					new_data_list = prev_data_list + selected_variant_CADD_data
					#print new_data_list
					annotated_dict[seq] = new_data_list

				#in order: chr, pos, ref, alt, conseq, consdetail, polyphen-cat, polyphen-val, CADD-raw, CADD-phred
				else:
					SNV_out +=1
					snv_keys_not_in_CADD + [snv_key]
					selected_variant_CADD_data = ['ref','ref','ref','ref','ref','ref','ref','ref','ref','ref',]
					annotated_dict[seq] = my_dict[seq]+selected_variant_CADD_data
		#add a second option for 1 bp deletions here
	print full_length_seqs, "full length seqs\t", single_nuc_subs, "SNVs"
	print SNV_in, "SNVs with CADD scores", SNV_out, "SNVs without CADD scores"
	print len(snv_keys_not_in_CADD), "snv_keys_not_in_CADD"
	#notes to self - why are so many SNVs not being found in the CADD score dictionary? -- possible off by 1 error/ indexing error in code?			
	return annotated_dict

def cigar_to_edits(cigar, seq, ref_seq):
	MDI_indexes = []
	for x in range(0,len(cigar)):
		if (cigar[x] == "M") or (cigar[x] == "D") or (cigar[x] == "I"):
			MDI_indexes.append(x)
	#now iterate over the MDI index list and do something different to call variants for each of the different options
	#these idexes will keep track of where we are in the read ('seq'), in the ref_seq, and the cigar
	ref_length = len(ref_seq)
	cigar_index = 0
	read_index = 0
	ref_index = 0
	#keeps matches listed
	mod_cigar = ''
	#only shows location of edits
	edit_string = ''
	for mdi_index in MDI_indexes:
		char = cigar[mdi_index]
		length = int(cigar[cigar_index:mdi_index])
		#cigar_index no longer useable for this chunk after this
		cigar_index = mdi_index+1
		#will call each SNV individually
		#X will be used to indicate a mismatch (an isolated X is a SNV - not specified here, though)
		# mod_cigar FORMAT: [bases][type][if insertion or mismatch, base id(s)]
		# edit_string FORMAT: comma-separated list of variants from WT:  [position]-[variant type]-[length (for I/D)]-[base id(s) (for I/X)]
		if char == 'M':
			match_length = 0
			#x will be relative to the reference
			#need an adjustment factor to get the appropriate index in the seq
			adj = ref_index-read_index
			for x in range(ref_index,ref_index+length):
				if seq[x-adj] == ref_seq[x]:
					match_length+=1
					if x == (ref_index+length-1):
						mod_cigar += (str(match_length)+'M')
						match_length = 0
				elif seq[x-adj] != ref_seq[x]:
					if match_length >= 1:
						mod_cigar += (str(match_length)+'M'+'1X'+seq[x-adj])
					elif match_length == 0:
						mod_cigar += ('1X'+seq[x-adj])
					edit_string += (str(x+1)+'-'+'X-'+seq[x-adj]+',')
					####### debug code:
					if x+1 > ref_length:
						print edit_string
						print ref_index, 'ref_index'
						print x, "x"
						print length, "length"
						print ref_index+x+1, "Ref_index+x+1"
						print ref_seq
						print seq
					match_length = 0
			#adjust the read_index and ref_index
			ref_index += length
			read_index += length
			
		elif char == 'D':
			mod_cigar += (str(length)+'D')
			edit_string += (str(ref_index+1)+'-D'+str(length)+',')
			#call position of the deleted bases
			#adjust the ref_index
			ref_index += length
		elif char == 'I':
			insertion_sequence = seq[read_index:read_index+length]
			mod_cigar += (str(length)+'I'+insertion_sequence)
			edit_string += (str(ref_index+1)+'-I'+str(length)+'-'+insertion_sequence+',')
			#adjust the read_index (ref_index stays the same)
			read_index += length
	if edit_string == '':
		edit_string = 'WT'
	else: #strip off the last comma
		edit_string = edit_string[:-1]
	#store the mod_cigar and the edit_string in the variant_dict
	#print ref_seq
	#print seq
	#print cigar
	#print mod_cigar
	#print edit_string

	#also return the number of edits, and the number of each type of edit

	return([mod_cigar, edit_string])

