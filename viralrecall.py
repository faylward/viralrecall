#!/usr/bin/env python
import sys, os, re, shlex, subprocess, pandas, numpy, itertools, argparse, time
from collections import defaultdict
from Bio import SeqIO
from operator import itemgetter
from itertools import islice
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

#pfam = "hmm/pfam.reduced.hmm"

# predict proteins from genome FNA file
def predict_proteins(genome_file, project, redo, batch):

	seqdict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
	if len(seqdict) < 1:
		raise (genome_file+" does not appear to be in FASTA format!")

	file_name = os.path.basename(genome_file)
	file_base = os.path.splitext(file_name)[0]
	if batch:
		path_base = os.path.splitext(project)[0]
		protein_file = os.path.join(path_base, file_base+".faa")
	else:
		base_name = os.path.basename(file_base)
		#protein_file = os.path.join(project, re.sub('.fna', '.faa', base_name))
		protein_file = os.path.join(project, project+".faa")	
	
	cmd = "prodigal -p meta -i "+ genome_file +" -a "+ protein_file
	#print(cmd)
	cmd2 = shlex.split(cmd)
	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
	return protein_file

	# get vog hmm descriptions
def get_annot(database):
	if database == "general":
		input = open("hmm/vog.annotations.tsv", "r")
		vdesc = defaultdict(lambda:"NA")
		for i in input.readlines():
			line = i.rstrip()
			tabs = line.split("\t")
			vog = tabs[0]
			desc = tabs[4]
			vdesc[vog] = desc
		return vdesc
	else:
		input = open("hmm/gvog_annotation.tsv", "r")
		vdesc = defaultdict(lambda:"NA")
		for i in input.readlines():
			line = i.rstrip("\n")
			tabs = line.split("\t")
			vog = tabs[0]
			desc = tabs[5]
			vdesc[vog] = desc
		return vdesc

# get accessions of HMM hits to exclude from bit score calculations
def get_accs(infile):
	acc2norm = defaultdict(lambda:float(1))
	handle = open(infile, "r")
	for i in handle.readlines():
		if i.startswith("ACC"):
			pass
		else:
			line = i.rstrip()
			tabs = line.split("\t")
			acc2norm[re.sub(".trim$", "", tabs[0])] = float(tabs[3])
			#acc_list.append(line)
	return(acc2norm)

# get SeqIO dictionary of input nucleic acid FASTA file
def get_fasta(genome_file):
	genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
	contig2length = {}
	for j in genome_dict:
		contig2length[j] = len(genome_dict[j])
	return genome_dict, contig2length

# get final prophage coordinates to output once desired flanking regions are accommodated
def get_finalcoords(contig_coords, subset_coords, flanking):
	newstart = int(subset_coords[0]) - int(flanking)
	if newstart < int(contig_coords[0]):
		newstart = int(0)

	newend = int(subset_coords[1]) + int(flanking)
	if newend > int(contig_coords[1]):
		newend = int(contig_coords[1])

	#print(contig_coords, subset_coords, flanking, newstart, newend)
	return((newstart, newend))	

# define sliding window function
def window(seq, n=25):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

# define function to acquire protein statistics from the faa file
def get_seqlist(input_file, project):
	prot2genome = {}
	record_dict = {}
	seqlist = {}
	strandlist = {}
	prot2start = {}
	prot2end = {}
	prot2contig = {}
	for i in SeqIO.parse(input_file, "fasta"):
		record_dict[i.id] = i
		seqlist[i.id] = len(i.seq)
		desc = i.description
		space = desc.split(" ")
		start = space[2]
		end = space[4]
		strand = space[6]
		strandlist[i.id] = strand
		prot2start[i.id] = start
		prot2end[i.id] = end
		prot2genome[i.id] = project

		protein = i.id
		protlist = protein.split("_")
		contig = "_".join(protlist[0:len(protlist)-1])
		prot2contig[protein] = contig

	return prot2genome, seqlist, strandlist, prot2start, prot2end, prot2contig, record_dict
		
# run HMMER3
def run_hmmer(input_file, db, suffix, cpus, redo, evalue):
	output_file = re.sub(".faa$", suffix, input_file)
	#print(output_file)
	if suffix == ".pfamout":
		cmd = "hmmsearch --cut_nc --cpu "+ cpus +" --tblout "+ output_file +" hmm/pfam.hmm "+ input_file
		#print(cmd)
	elif suffix == ".vogout":
		if db == "general":
			vogdb = "hmm/vogdb.hmm"
			cmd = "hmmsearch --cpu "+ cpus +" --tblout "+ output_file +" "+ vogdb +" "+ input_file
		elif db == "GVOG":
			vogdb = "hmm/gvog.hmm"
			cmd = "hmmsearch --cpu "+ cpus +" --tblout "+ output_file +" "+ vogdb +" "+ input_file	
			#print(cmd)
		elif db == "marker":
			vogdb = "hmm/NCLDV_markers.hmm"
			cmd = "hmmsearch --cpu "+ cpus +" --tblout "+ output_file +" "+ vogdb +" "+ input_file	
#			#print(cmd)

	#print(db)
	cmd2 = shlex.split(cmd)

	if not redo and not db == "marker":
		#pass
		#print(cmd)
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
		#os.remove("out.txt")
	return output_file

# define function for parsing HMMER3 output
def parse_hmmout(hmmout, evalue):

	input = open(re.sub("\/$", "", hmmout), "r")
	hit_dict = defaultdict(lambda:"NA")
	bit_dict = defaultdict(float)
	hit2pfam = {}
	gvog2norm = get_accs("acc/gvog_summary.tsv")
	pfam2norm = get_accs("acc/pfam_summary.tsv")

	for i in input.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			protein = tabs[0]
			hit = re.sub(".trim$", "", tabs[2])
			pfamhit = tabs[3]
			#print pfamhit
			eval = float(tabs[4])
			#print pfamhit, eval

			if float(tabs[5]) > 0:
				if "vogout" in hmmout:
					score = (math.sqrt(float(tabs[5]))) * gvog2norm[hit]
					#score = float(tabs[5]) * gvog2norm[hit]
					#print (protein, tabs[5], math.sqrt(float(tabs[5])), gvog2norm[hit], score, hit)
				else:
					score = (math.sqrt(float(tabs[5]))) * pfam2norm[pfamhit]
					#score = float(tabs[5]) * pfam2norm[pfamhit]
					#print(protein, pfamhit, tabs[5], pfam2norm[pfamhit], score)
			else:
				score = 0

			if "vogout" in hmmout:
				if score > bit_dict[protein] and eval <= float(evalue):
					bit_dict[protein] = score
					hit_dict[protein] = hit
					hit2pfam[protein] = pfamhit
				else:
					pass

			else:
				if score > bit_dict[protein]:
					bit_dict[protein] = score
					hit_dict[protein] = hit
					hit2pfam[protein] = pfamhit
				else:
					pass

	return hit_dict, bit_dict

# run HMMER3 on marker genes only
def marker_hmmer(input_file, cpus, redo):
	output_file = re.sub(".faa", ".markerout", input_file)
	cmd = "hmmsearch --cpu "+ cpus +" --tblout "+ output_file +" hmm/NCLDV_markers.hmm "+ input_file	
	cmd2 = shlex.split(cmd)

	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
		#os.remove("out.txt")
	return output_file

def parse_markers(hmmout):
	input = open(hmmout, "r")
	score_dict = {"A32":float(80), "D5":float(80), "SFII":float(100), "mcp":float(80), "mRNAc":float(80), "PolB":float(200), "RNAPL":float(200), "RNAPS":float(200), "RNR":float(80), "VLTF3":float(80)}
	contig2hits = defaultdict(list)
	protein2hits = defaultdict(lambda:"-")

	for i in input.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			protein = tabs[0]
			hit = tabs[2]
			score = float(tabs[5])
			contig = re.sub("_\d+$", "", protein)
			
			if score > score_dict[hit]:
				#print(protein, contig, hit, score)
				contig2hits[contig].append(hit +":"+ str(score))
				protein2hits[protein] = hit +":"+ str(score)
				
	contig2final = defaultdict(lambda:"-")
	for n in contig2hits:
		outstr = ",".join(contig2hits[n])
		contig2final[n] = outstr
	#print(contig2final)
	return contig2final, protein2hits
	
	
# define function to get continuous genome coordinates in case of multiple contigs/plasmids
def cumsum2(list1, contig2length, prot2start):
	baseval = 0
	cumsum = []
	already_done = []
	for index,i in enumerate(list1):
		start = float(prot2start[i])
		contig = re.sub("_\d+$", "", i)
		if len(already_done) < 1:
			cumsum.append(start)
			
			if contig in already_done:
				pass
			else:
				already_done.append(contig)
				
		else:
			if contig in already_done:
				pass
			else:
				already_done.append(contig)
			
			#print(already_done)
			baseval = 0
			for j in already_done[0:len(already_done)-1]:
				baseval += contig2length[j]
		
			new = baseval + start
			cumsum.append(new)
			#print(new)
	return cumsum
		#print(index, i)
	#	if index == 0:
	#		cumsum.append(start)
			
	#	elif start < float(prot2start[list1[index-1]]):
	#		contig = re.sub("_\d+$", "", i)
			#print(contig, contig2length[contig])
			#baseval = cumsum[index-1]
	#		baseval = contig2length[contig] + baseval
	#		new = start + baseval
	#		cumsum.append(new)
	#		print(i, prot2start[i], new, baseval)
			#print(i, baseval)
	#	else:
	#		new = start + baseval
	#		cumsum.append(new)
	#return cumsum

# get genomic regions that look like phage
def get_regions(list1):
	index_list = []
	for index,value in enumerate(list1):
		if index == 0:
			pass
		#elif value == 0:
		#	pass
		elif value >= 0 and list1[index-1] >= 0:
			newval = index-1
			if newval in index_list:
				pass
			else:
				index_list.append(index-1)
			index_list.append(index)
	return index_list

# main function that runs the program
def run_program(input, project, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel):

	# create output directories
	foldername = os.path.splitext(project)[0]
	if os.path.isdir(foldername):
		#raise Exception(project+' already exists')
		#print("\n*******************************************************************************\nOverwriting existing project folder! Waiting 5 seconds so you have time to exit\n*******************************************************************************")
		#time.sleep(5)
		pass
	else:
		os.mkdir(foldername)

	# remove previous files before re-calculating results
	if redo:
		for files in os.listdir(foldername):
			if files.endswith(".tsv") or "_viral_region_" in files:
				os.remove(os.path.join(foldername, files))

	if batch:
		relpath = os.path.split(project)[1]
		relpathbase = os.path.splitext(relpath)[0]
		base = os.path.splitext(project)[0]
		#print(project, base, relpath, relpathbase)
		
	# predict proteins, run HMMER3 searches, and parse outputs
	#print(input, project, redo, batch)
	protein_file = predict_proteins(input, project, redo, batch)

	if database == "marker":
		#vog_out  = run_hmmer(protein_file, database, ".vogout", cpus, redo, evalue)
		#pfam_out = run_hmmer(protein_file, "hmm/pfam.hmm", ".pfamout", cpus, redo, evalue)		

		pfam_hit = defaultdict(lambda:"no-search-performed")
		pfam_bit = defaultdict(lambda:float(0))

		vog_hit = defaultdict(lambda:"no-search-performed")
		vog_bit = defaultdict(lambda:float(0))
		
	else:
	
		if not os.path.exists("hmm/gvog.hmm") or not os.path.exists("hmm/vogdb.hmm"):
			raise("Can't seem to find the hmm databases in the hmm/ directory. Please see GitHub for download instructions")
		
		vog_out  = run_hmmer(protein_file, database, ".vogout", cpus, redo, evalue)
		pfam_out = run_hmmer(protein_file, "hmm/pfam.hmm", ".pfamout", cpus, redo, evalue)

		vog_hit, vog_bit = parse_hmmout(vog_out, evalue)
		pfam_hit, pfam_bit = parse_hmmout(pfam_out, evalue)

		
	marker_out = marker_hmmer(protein_file, cpus, redo)
	contig2markers, protein2markers = parse_markers(marker_out)
		
	vdesc = get_annot(database)
	vog_annot = {}
	for i in vog_hit:
		vog_annot[i] = vdesc[vog_hit[i]]
		


	# get protein features from Prodigal FASTA headers
	prot2genome, seqs, strands, prot2start, prot2end, prot2contig, record_dict = get_seqlist(protein_file, os.path.basename(project))

	# get dictionary of nucleic acid sequences
	genome_dict, contig2length = get_fasta(input)
	
	# set up Pandas DataFrame of the full genome annotation
	names = ['genome', 'replicon', 'vog', 'virbit', 'vdesc', 'pfam', 'pfambit', 'protlength', 'strand', 'start', 'end']
	df = pandas.DataFrame()
	for index, i in enumerate([prot2genome, prot2contig, vog_hit, vog_bit, vog_annot, pfam_hit, pfam_bit, seqs, strands, prot2start, prot2end]):
		s1 = pandas.DataFrame(pandas.Series(i, name = names[index]))
		df = pandas.concat([df, s1], axis=1, sort=True)

	#print(df)
	# fill all NA values with 0, ensuring that proteins with no hit to either Pfam or VOG are counted as having a bit score of 0. The "prophage score" is then calculated as the difference between the Pfam and VOG scores.
	df.fillna(float(0), inplace=True, axis=1)
	df["start"] = pandas.to_numeric(df['start'])
	df["end"] = pandas.to_numeric(df['end'])
	df = df.sort_values(by=['replicon', 'start'])
	df["score"] = df["virbit"] - df["pfambit"]

	# now for the plot we need a contantly increasing axis so we don't plot contigs/plasmids over each other
	starts = df["start"].tolist()
	starts = [float(i) for i in starts]
	cumsum = []
	#df["cumsum"] = cumsum2(df["start"], contig2length, prot2start)
	df["cumsum"] = cumsum2(df.index, contig2length, prot2start)
	
	# for each replicon we need to go through and calculate a rolling mean of the prophage scores. 
	contigs = sorted(set(prot2contig.values()))
	contig_bounds = []
	df2 = pandas.DataFrame()
	for index, contig in enumerate(contigs):
		# we need to calculate the rolling mean separately for each replicon so we don't get overlap between non-contiguous sequences
		subset = pandas.DataFrame(df.loc[df['replicon'] == contig])
		subset["rolling"] = subset["score"].rolling(window, min_periods=3, center=True).mean()
		df2 = pandas.concat([df2, subset])

		ends = subset["end"].tolist()
		ends = [float(i) for i in ends]
		contig_end = max(ends)

		#print(contig_end, index)
		# the contig bounds are used later for plotting, so we know where replicons end
		if len(contig_bounds) > 0:
			new_end = float(contig_end) + float(contig_bounds[index-1])
		else:
			new_end = contig_end
		contig_bounds.append(new_end)

	#df2 = df2.sort_values(by=["cumsum"])
	df2.fillna(0, inplace=True, axis=1)

	reps = set(df2["replicon"].tolist())

	# initialize summary dataframe that we will append to as we find viral regions
	summary = pandas.DataFrame()
	tally = 0

	# proceed here if all you want is the contig-level stats (no info about viral regions)
	if contiglevel:
		for rep in reps:
			df3 = df2.loc[df2['replicon'] == rep]
			tally +=1

			#indices = list(map(itemgetter(1), group))
			#subset = df3.ix[indices]
			minval = min(df3["start"])
			maxval = max(df3["end"])

			vogacc = df3["vog"].tolist()
			num_prot = len(df3["protlength"].tolist())

			score = numpy.mean(df3["score"])
			voghits = len([i for i in df3["virbit"].tolist() if i > 0])
			pfamhits = len([i for i in df3["pfambit"].tolist() if i > 0])
			length = int(float(maxval) - float(minval))
			replicon = rep

			record = genome_dict[replicon]
			contig_length = len(record.seq)

			markerhits = contig2markers[rep]
			
			data = pandas.Series([replicon, contig_length, score, num_prot, voghits, pfamhits, markerhits], name=project)
			summary = summary.append(data)

		# output summary files
		if batch:
			df2["vog"] = df2["vog"].replace(0, "no_hit")
			df2["vdesc"] = df2["vdesc"].replace(0, "no_hit")
			df2["pfam"] = df2["pfam"].replace(0, "no_hit")
			df2.to_csv(os.path.join(base, relpathbase+".full_annot.tsv"), sep="\t", index_label="protein_ids")

			if summary.shape[1] > 0:

				summary.columns = ['replicon', 'contig_length', 'score', 'num_ORFs', 'num_viralhits', 'num_pfamhits', 'markerhits']
				#base = os.path.basename(project)
				summary_file.write(base +"\t"+ str(summary.shape[0]) +"\n")
				summary.to_csv(os.path.join(base, relpathbase+".summary.tsv"), sep="\t", index_label="viral_regions")

			else:
				#base = os.path.basename(project)
				summary_file.write(base +"\t0\n")

		else:
			if summary.shape[1] > 0:
				summary.columns = ['replicon', 'contig_length', 'score', 'num_ORFs', 'num_viralhits', 'num_pfamhits', 'markerhits']
				summary.to_csv(os.path.join(project, project+".summary.tsv"), sep="\t", index_label="project")
				
				if database == "marker":
					df2["vog"] = df2["vog"].replace(0, "no-search-performed")
					df2["vdesc"] = df2["vdesc"].replace(0, "no-search-performed")
					df2["pfam"] = df2["pfam"].replace(0, "no-search-performed")
				else:
					df2["vog"] = df2["vog"].replace(0, "no_hit")
					df2["vdesc"] = df2["vdesc"].replace(0, "no_hit")
					df2["pfam"] = df2["pfam"].replace(0, "no_hit")
				
				df2.to_csv(os.path.join(project, project+".full_annot.tsv"), sep="\t", index_label="protein_ids")

	# otherwise proceed with regular viralrecall to identify virus-like regions 
	else:
		for rep in reps:

			df3 = df2.loc[df2['replicon'] == rep]
			#print(rep, df3.shape)
			# now let's get the regions of the entire genome file that have a net positive prophage signal
			reg = get_regions(df3["rolling"].tolist())
			reg = [int(i) for i in reg]

			# now let's subset the genome to get only the prophage regions, and output that so we can look at it later if we want
			#subset = df3.ix[reg]
			subset = df3.iloc[reg]
			if batch:
				#print(os.path.join(project, base+".vregion_annot.tsv"), project, base)
				subset.to_csv(os.path.join(base, relpathbase+".vregion_annot.tsv"), sep='\t', index_label="protein_ids")
			else:
				subset.to_csv(os.path.join(project, project+".vregion_annot.tsv"), sep='\t', index_label="protein_ids")

			# now let's get a summary of each prophage region, and output that
			for key, group in itertools.groupby(enumerate(reg), key=lambda ix:ix[0]-ix[1]):
				indices = list(map(itemgetter(1), group))
				#subset = df3.ix[indices]
				subset = df3.iloc[indices]
				minval = min(subset["start"])
				maxval = max(subset["end"])
				#print(minval, maxval)
				#print(key, map(itemgetter(1), group), group, indices, [replicons[k] for k in indices])
				#print([replicons[k] for k in indices])

				vogacc = subset["vog"].tolist()

				score = numpy.mean(subset["score"])
				voghits = len([i for i in subset["virbit"].tolist() if i > 0])
				length = int(float(maxval) - float(minval))
				replicon = subset['replicon'].tolist()[0]
				markers = [protein2markers[j] for j in list(subset.index) if not protein2markers[j] == "-"]
				markerlist = ",".join(markers)
				#print(minval, maxval, replicon, score, length)
		
				# Let's filter the putatige prophage by the parameters used in the input. 
				if length >= phagesize and score > minscore and voghits >= minhit:
					#print(minval, maxval, replicon, score, length)
					tally +=1

					record = genome_dict[replicon]
					contig_length = len(record.seq)

					data = pandas.Series([replicon, minval, maxval, length, contig_length, score, voghits, len(indices), markerlist], name="viral_region_"+str(tally))
					summary = summary.append(data)
					#print(summary)

					# now let's output the proteins and nucleic acid sequence of the putative prophage
					if batch:
						#base = os.path.basename(project)
						protein_file = os.path.join(base, relpathbase+"_viral_region_"+str(tally)+".faa")
					else:
						protein_file = os.path.join(project, project+"_viral_region_"+str(tally)+".faa")

					proteins = list(subset.index)
					records = [record_dict[record] for record in record_dict.keys() if record in proteins]
					SeqIO.write(records, protein_file, "fasta")

					if batch:
						#base = os.path.basename(project)
						nucl_file = open(os.path.join(base, relpathbase+"_viral_region_"+str(tally)+".fna"), "w")
					else:
						nucl_file = open(os.path.join(project, project+"_viral_region_"+str(tally)+".fna"), "w")

					record = genome_dict[replicon]
					seq = record.seq
			
					newcoords = get_finalcoords((0, len(record.seq)), (minval, maxval), flanking)
					prophage_region = seq[newcoords[0]:newcoords[1]]
					nucl_file.write(">"+ project+"_viral_region_"+str(tally) +" "+ record.id +"\n"+ str(prophage_region))
			
			# if we find any prophage let's output a summary file
		if batch:
			df2["vog"] = df2["vog"].replace(0, "no_hit")
			df2["vdesc"] = df2["vdesc"].replace(0, "no_hit")
			df2["pfam"] = df2["pfam"].replace(0, "no_hit")
			df2.to_csv(os.path.join(base, relpathbase+".full_annot.tsv"), sep="\t", index_label="protein_ids")

			if summary.shape[1] > 0:

				summary.columns = ['replicon', 'start_coord', 'end_coord', 'vregion_length', 'contig_length', 'score', 'num_viralhits', 'num_ORFs', 'markers']
				#base = os.path.basename(project)
				summary_file.write(base +"\t"+ str(summary.shape[0]) +"\n")
				summary.to_csv(os.path.join(base, relpathbase+".summary.tsv"), sep="\t", index_label="viral_regions")

			else:
				#base = os.path.basename(project)
				summary_file.write(base +"\t0\n")

		else:
			if summary.shape[1] > 0:
				summary.columns = ['replicon', 'start_coord', 'end_coord', 'vregion_length', 'contig_length', 'score', 'num_viralhits', 'num_ORFs', 'markers']
				summary.to_csv(os.path.join(project, project+".summary.tsv"), sep="\t", index_label="viral_regions")
			else:
				summary_out = open(os.path.join(project, project+".summary.tsv"), "w")
				summary_out.write('replicon\tstart_coord\tend_coord\tvregion_length\tcontig_length\tscore\tnum_viralhits\tnum_ORFs\tmarkers\n')


			if database == "marker":
				df2["vog"] = df2["vog"].replace(0, "no-search-performed")
				df2["vdesc"] = df2["vdesc"].replace(0, "no-search-performed")
				df2["pfam"] = df2["pfam"].replace(0, "no-search-performed")
			else:
				df2["vog"] = df2["vog"].replace(0, "no_hit")
				df2["vdesc"] = df2["vdesc"].replace(0, "no_hit")
				df2["pfam"] = df2["pfam"].replace(0, "no_hit")

			df2.to_csv(os.path.join(project, project+".full_annot.tsv"), sep="\t", index_label="protein_ids")


	#######################################################
	################# Print figure ########################
	if (plotflag):
		if database == "marker":
			print("Plotting options disabled when db = marker!")
		else:	
			f = plt.figure(figsize=(15,4))

			maxbound = round(max(df2["cumsum"])/1000000) * 1000000
			maxbound_label = int(maxbound / 1000000)
			
			bound_labels = list(range(int(0), maxbound_label+1))
			bounds = [item*1000000 for item in bound_labels]
			
			val = numpy.nanmax(df2["rolling"])
			minval = numpy.nanmin(df2["rolling"])
			#print(bound_labels, bounds, contig_bounds, val, minval)
			if len(contig_bounds) > 1:
				#print(contig_bounds, bound_labels)
				plt.vlines(contig_bounds, 0, val, colors="grey", zorder=20, linestyles = "dotted") #, linewidth=0.9)
				
			plt.plot(df2["cumsum"], df2["rolling"], color="grey")
			#plt.xticks([])
			plt.xlabel("Genome Position")
			plt.ylabel("Score")

			plt.fill_between(df2["cumsum"], df2["rolling"], where=(df2["rolling"] >= 0), facecolor="dodgerblue", alpha=0.8, interpolate=True)
			plt.fill_between(df2["cumsum"], df2["rolling"], where=(df2["rolling"] < 0), facecolor="firebrick", alpha=0.8, interpolate=True)
			#plt.xticks(bounds, bound_labels)
			plt.ylim(minval, numpy.nanmax(df2["rolling"]))
			if batch:
				f.savefig(os.path.join(base, relpathbase+".pdf"), bbox_inches='tight')
			else:
				f.savefig(os.path.join(project, project+".pdf"), bbox_inches='tight')
			plt.close()
	#######################################################

########################################################################
##### use argparse to run through the command line options given #######
########################################################################
def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="ViralRecall v. 2.0: A flexible command-line tool for predicting NCLDV-like regions in genomic data \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=True, help='Input FASTA file (ending in .fna)')
	args_parser.add_argument('-p', '--project', required=True, help='project name for outputs')
	args_parser.add_argument('-db', '--database', required=False, default="GVOG", help='Viral HMM database to use. Options are "general" for the general VOG db, "GVOG" for the GVOG db, and "marker" for searching only a set of 10 conserved NCLDV markers (good for screening large datasets). See README for details')
	args_parser.add_argument('-w', '--window', required=False, default=int(15), help='sliding window size to use for detecting viral regions (default=15)')
	args_parser.add_argument('-m', '--minsize', required=False, default=int(10), help='minimum length of viral regions to report, in kilobases (default=10)')
	args_parser.add_argument('-s', '--minscore', required=False, default=int(1), help='minimum score of viral regions to report, with higher values indicating higher confidence (default=1)')
	args_parser.add_argument('-g', '--minhit', required=False, default=int(4), help='minimum number of viral hits that each viral region must have to be reported (default=4)')
	args_parser.add_argument('-e', '--evalue', required=False, default=str(1e-10), help='e-value that is passed to HMMER3 for the VOG hmmsearch (default=1e-10)')
	args_parser.add_argument('-fl', '--flanking', required=False, default=int(0), help='length of flanking regions upstream and downstream of the viral region to output in the final .fna files (default=0)')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	args_parser.add_argument('-b', '--batch', type=bool, default=False, const=True, nargs='?', help='Batch mode: implies the input is a folder of .fna files that each will be run iteratively')
	args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	args_parser.add_argument('-c', '--contiglevel', type=bool, default=False, const=True, nargs='?', help='calculate contig/replicon level statistics instead of looking at viral regions (good for screening contigs)')
	args_parser.add_argument('-f', '--figplot', type=bool, default=False, const=True, nargs='?', help='Specify this flag if you would like a plot of the viral-like regions with the output')
	args_parser.add_argument('-v', '--version', action='version', version='ViralRecall v. 2.1')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	input = args_parser.input
	project = args_parser.project
	database = args_parser.database
	window = int(args_parser.window)
	phagesize = int(args_parser.minsize)*1000
	minscore = int(args_parser.minscore)
	minhit = int(args_parser.minhit)
	evalue = str(args_parser.evalue)
	cpus = args_parser.cpus
	plotflag = args_parser.figplot
	redo = args_parser.redo
	contiglevel = args_parser.contiglevel
	flanking = args_parser.flanking
	batch = args_parser.batch

	project = project.rstrip("/")
	if batch:
		if os.path.isdir(project):
			pass
		else:
			os.mkdir(project)
		
		summary_file = open(os.path.join(project, "batch_summary.txt"), "w")
		summary_file.write("genome\tcontigs_tested\n")

		if os.path.isdir(project):
			pass
		else:
			os.mkdir(project)
			
		file_list = os.listdir(input)
		for i in file_list:
			#if i.endswith(".fna"):
			#name = re.sub(".fna", "", i)
			newproject = os.path.join(project, i)
			#newproject = os.path.splitext(newproject)[0]
			newinput = os.path.join(input, i)
			print("Running viralrecall on "+ i + " and output will be deposited in "+ newproject)
			run_program(newinput, newproject, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel)
	else:
		summary_file = 1
		run_program(input, project, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

# end




