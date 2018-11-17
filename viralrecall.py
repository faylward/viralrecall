import sys, os, re, shlex, subprocess, pandas, numpy, itertools, argparse, time
from collections import defaultdict
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from operator import itemgetter
from itertools import islice

# locations for the HMM database files to be used
vogdb = "hmm/vogdb.hmm"
pfam = "hmm/pfam.reduced.hmm"

# predict proteins from genome FNA file
def predict_proteins(genome_file, project, redo):
	protein_file = os.path.join(project, re.sub(".fna", ".faa", genome_file))
	cmd = "prodigal -i "+ genome_file +" -a "+ protein_file
	cmd2 = shlex.split(cmd)
	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
	return protein_file

# get SeqIO dictionary of input nucleic acid FASTA file
def get_fasta(genome_file):
	genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
	return genome_dict		

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
def get_seqlist(input_file):
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

		protein = i.id
		protlist = protein.split("_")
		contig = "_".join(protlist[0:len(protlist)-1])
		prot2contig[protein] = contig

	return seqlist, strandlist, prot2start, prot2end, prot2contig, record_dict
		
# run HMMER3
def run_hmmer(input_file, db, suffix, cpus, redo):
	output_file = re.sub(".faa", suffix, input_file)
	if suffix == ".pfamout":
		cmd = "hmmsearch --cpu "+ cpus +" --cut_nc --tblout "+ output_file +" "+ db +" "+ input_file
	else:
		cmd = "hmmsearch --cpu "+ cpus +" -E 1e-10 --tblout "+ output_file +" "+ db +" "+ input_file	
	#print(cmd)
	cmd2 = shlex.split(cmd)
	if not redo:
		subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))
	return output_file

# define function for parsing HMMER3 output
def parse_hmmout(hmmout):
	input = open(hmmout, "r")
	hit_dict = defaultdict(lambda:"NA")
	bit_dict = defaultdict(float)

	for i in input.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			protein = tabs[0]
			hit = tabs[2]
			eval = float(tabs[4])
			score = float(tabs[5])

			if score > bit_dict[protein]:
				bit_dict[protein] = score
				hit_dict[protein] = hit
			else:
				pass
	return hit_dict, bit_dict

# define function to get continuous genome coordinates in case of multiple contigs/plasmids
def cumsum2(list1):
	baseval = 0
	cumsum = []
	for index,i in enumerate(list1):
		if index==0:
			cumsum.append(i)
		elif i < list1[index-1]:
			baseval = cumsum[index-1]
			new = i + baseval
			cumsum.append(new)
			#print(i, baseval)
		else:
			new = i + baseval
			cumsum.append(new)
	return cumsum

# get genomic regions that look like phage
def get_regions(list1):
	index_list = []
	for index,value in enumerate(list1):
		if index == 0:
			pass
		elif value == 0:
			pass
		elif value > 0 and list1[index-1] > 0:
			index_list.append(index-1)
	return index_list

# main function that runs the program
def run_program(input, project, window, phagesize, minscore, minvog, cpus, plotflag, redo, flanking):

	# create output directories
	if os.path.isdir(project):
		#raise Exception(project+' already exists')
		#print("\n*******************************************************************************\nOverwriting existing project folder! Waiting 5 seconds so you have time to exit\n*******************************************************************************")
		#time.sleep(5)
		pass
	else:
		os.mkdir(project)

	# predict proteins, run HMMER3 searches, and parse outputs
	protein_file = predict_proteins(input, project, redo)
	vog_out  = run_hmmer(protein_file, vogdb, ".vogout", cpus, redo)
	pfam_out = run_hmmer(protein_file, pfam, ".pfamout", cpus, redo)

	vog_hit, vog_bit = parse_hmmout(vog_out)
	pfam_hit, pfam_bit = parse_hmmout(pfam_out)

	# get protein features from Prodigal FASTA headers
	seqs, strands, prot2start, prot2end, prot2contig, record_dict = get_seqlist(protein_file)

	# get dictionary of nucleic acid sequences
	genome_dict = get_fasta(input)
	
	# set up Pandas DataFrame of the full genome annotation
	names = ['replicon', 'vog', 'vogbit', 'pfam', 'pfambit', 'protlength', 'strand', 'start', 'end']
	df = pandas.DataFrame()
	for index, i in enumerate([prot2contig, vog_hit, vog_bit, pfam_hit, pfam_bit, seqs, strands, prot2start, prot2end]):
		s1 = pandas.DataFrame(pandas.Series(i, name = names[index]))
		df = pandas.concat([df, s1], axis=1, sort=True)

	# fill all NA values with 0, ensuring that proteins with no hit to either Pfam or VOG are counted as having a bit score of 0. The "prophage score" is then calculated as the difference between the Pfam and VOG scores.
	df.fillna(float(0), inplace=True, axis=1)
	df["start"] = pandas.to_numeric(df['start'])
	df = df.sort_values(by=['replicon', 'start'])
	df["score"] = df["vogbit"] - df["pfambit"]

	# now for the plot we need a contantly increasing axis so we don't plot contigs/plasmids over each other
	starts = df["start"].tolist()
	starts = [float(i) for i in starts]
	cumsum = []
	df["cumsum"] = cumsum2(df["start"])

	# for each replicon we need to go through and calculate a rolling mean of the prophage scores. 
	contigs = sorted(set(prot2contig.values()))
	contig_bounds = []
	df2 = pandas.DataFrame()
	for index, contig in enumerate(contigs):
		# we need to calculate the rolling mean separately for each replicon so we don't get overlap between non-contiguous sequences
		subset = pandas.DataFrame(df.loc[df['replicon'] == contig])
		subset["rolling"] = subset["score"].rolling(window, min_periods=3).mean()
		df2 = pandas.concat([df2, subset])

		ends = subset["end"].tolist()
		ends = [float(i) for i in ends]
		contig_end = max(ends)

		# the contig bounds are used later for plotting, so we know where replicons end
		if len(contig_bounds) > 0:
			new_end = float(contig_end) + float(contig_bounds[index-1])
		else:
			new_end = contig_end
		contig_bounds.append(new_end)

	df2 = df2.sort_values(by=["cumsum"])
#	df2.fillna(0, inplace=True, axis=1)

	# now let's get the regions of the entire genome file that have a net positive prophage signal
	reg = get_regions(df2["rolling"].tolist())
	reg = [int(i) for i in reg]

	# now let's subset the genome to get only the prophage regions, and output that so we can look at it later if we want
	subset = df2.ix[reg]
	subset.to_csv(os.path.join(project, project+".prophage_annot.tsv"), sep='\t', index_label="protein_ids")

	# now let's get a summary of each prophage region, and output that
	tally = 0
	summary = pandas.DataFrame()
	for key, group in itertools.groupby(enumerate(reg), key=lambda ix:ix[0]-ix[1]):
		indices = list(map(itemgetter(1), group))
		#if len(indices) >= phagesize:
		subset = df.ix[indices]
		minval = min(subset["start"])
		maxval = max(subset["start"])
		
		score = numpy.mean(subset["score"])
		voghits = len([i for i in subset["vogbit"].tolist() if i > 0])
		length = int(maxval - minval)
		replicon = subset['replicon'].tolist()[0]

		# Let's filter the putatige prophage by the parameters used in the input. 
		if length >= phagesize and score > minscore and voghits > minvog:
			tally +=1
			data = pandas.Series([replicon, minval, maxval, length, score, voghits, len(indices)], name="prophage_Region_"+str(tally))
			summary = summary.append(data)

			# now let's output the proteins and nucleic acid sequence of the putative prophage
			protein_file = os.path.join(project, project+"_prophage_region_"+str(tally)+".faa")
			proteins = list(subset.index)
			records = [record_dict[record] for record in record_dict.keys() if record in proteins]
			SeqIO.write(records, protein_file, "fasta")

			nucl_file = open(os.path.join(project, project+"_prophage_region_"+str(tally)+".fna"), "w")
			record = genome_dict[replicon]
			seq = record.seq
			
			newcoords = get_finalcoords((0, len(record.seq)), (minval, maxval), flanking)
			prophage_region = seq[newcoords[0]:newcoords[1]]
			nucl_file.write(">"+ project+"_prophage_region_"+str(tally) +" "+ record.id +"\n"+ str(prophage_region))
			
	# if we find any prophage let's output a summary file
	if summary.shape[1] > 0:
		summary.columns = ['replicon', 'start_coord', 'end_coord', 'prophage_length', 'score', 'num_voghits', 'num_ORFs']
		summary.to_csv(os.path.join(project, project+".summary.tsv"), sep="\t", index_label="prophage_regions")
		df2.to_csv(os.path.join(project, project+".full_annot.tsv"), sep="\t", index_label="protein_ids")

	#######################################################
	################# Print figure ########################
	if (plotflag):
		f = plt.figure(figsize=(15,4))

		maxbound = round(max(df2["cumsum"])/1000000) * 1000000
		maxbound_label = int(maxbound / 1000000)

		bound_labels = list(range(int(0), maxbound_label+1))
		bounds = [item*1000000 for item in bound_labels]

		val = max(df2["rolling"])/10
		plt.vlines(contig_bounds, 0, val, colors="red", zorder=20)
		plt.plot(df2["cumsum"], df2["rolling"])
		plt.xticks([])
		plt.xlabel("Genome Position (Mbp)")
		plt.ylabel("Score")

		plt.fill_between(df2["cumsum"], df2["rolling"])
		plt.xticks(bounds, bound_labels)
		plt.ylim(0, numpy.nanmax(df2["rolling"]))
		f.savefig(os.path.join(project, project+".pdf"), bbox_inches='tight')
	#######################################################

########################################################################
##### use argparse to run through the command line options given #######
########################################################################
def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="ViralRecall: Predicting prophage-like regions in prokaryotic genomes \nFrank O. Aylward, Assistant Professor, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='!!!\nIf you use this tool in a publication please do not forget to cite:\nProdigal (DOI 10.1186/1471-2105-11-119)\nHMMER3 (DOI 10.1371/journal.pcbi.1002195)\n!!!')
	args_parser.add_argument('-i', '--input', required=True, help='Input FASTA file (ending in .fna)')
	args_parser.add_argument('-p', '--project', required=True, help='project name for outputs')
	args_parser.add_argument('-w', '--window', required=False, default=int(15), help='sliding window size to use for detecting prophage regions (default=15)')
	args_parser.add_argument('-m', '--minsize', required=False, default=int(10), help='minimum length of prophage to report, in kilobases (default=10)')
	args_parser.add_argument('-s', '--minscore', required=False, default=int(20), help='minimum score of prophage to report, with higher values indicating higher confidence (default=20)')
	args_parser.add_argument('-v', '--minvog', required=False, default=int(4), help='minimum number of VOG hits that each prophage must have to be reported (default=4)')
	args_parser.add_argument('-fl', '--flanking', required=False, default=int(0), help='length of flanking regions upstream and downstream of the prophage to output in the final .fna files (default=0)')
	args_parser.add_argument('-t', '--cpus', required=False, default=1, help='number of cpus to use for the HMMER3 search')
#	args_parser.add_argument('-b', '--batch', type=bool, default=False, const=True, nargs='?', help='Batch mode: implies the input is a folder of .fna files that each will be run iteratively')
	args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	args_parser.add_argument('-f', '--figplot', type=bool, default=False, const=True, nargs='?', help='Specify this flag if you would like a plot of the prophage-like regions with the output')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	input = args_parser.input
	project = args_parser.project
	window = int(args_parser.window)
	phagesize = int(args_parser.minsize)*1000
	minscore = int(args_parser.minscore)
	minvog = int(args_parser.minvog)
	cpus = args_parser.cpus
	plotflag = args_parser.figplot
	redo = args_parser.redo
	flanking = args_parser.flanking
#	batch = args_parser.batch

#	if batch:
#		os.mkdir(project)
#		file_list = os.listdir(input)
#		for i in file_list:
#			if i.endswith(".fna"):
#				name = re.sub(".fna", "", i)
#				project = os.path.join(project, name)
#				input = os.path.join(input, i)
#				print(i, input, project, window, phagesize)
#				run_program(input, project, window, phagesize, minscore, cpus, plotflag, redo, batch)
#	else:
	run_program(input, project, window, phagesize, minscore, minvog, cpus, plotflag, redo, flanking)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

# end



