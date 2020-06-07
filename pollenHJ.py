#!/usr/bin/env python3

import sys
import argparse
import os
import textwrap
import pysam
import time, threading

def version():
	v = """
    ##########################################################################################
 
                                 pollen genetic mapping pipeline
 
    ------------------------------------------------------------------------------------------
    
                                       Version : 0.0.1.8
                                        Author : Yumin Huang

          Update:
             1. bam file NM tag added

    ##########################################################################################
	"""
	return v

def warning(*objs):
	print("WARNING: ", *objs, end='\n', file=sys.stderr)
	sys.exit()

def get_parser():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=textwrap.dedent(version())
	)

	group1 = parser.add_argument_group('Main', 'main arguments')
	group1.add_argument('-i','--input', help='bam input file, must be indexed', type=str)
	group1.add_argument('-s','--snps', help='SNPs info betweent parents', type=str)
	group1.add_argument('-o', '--output', help='output name', type=str)
	group1.add_argument('-f', '--fai', help='genome index file', type=str)
	group1.add_argument('-q', '--mapq', help='minimum MAPQ value (optional)', type=int, default=0)
	group1.add_argument('-n', '--NM', help='minimum mismatch number (default=10)', type=int, default=10)
	group1.add_argument('-Q', '--qs', help='minimum base qualities (optional)', type=int, default=0)
	group1.add_argument('-w', '--window', help='window size for each bin (default=2000000)', type=int, default=2000000)
	group1.add_argument('-e', '--step', help='window step for bins (default=100000)', type=int, default=100000)
	group1.add_argument('--insertion', help='insertion fragment size of sequencing', type=int, default=400)
	group1.add_argument('--tp', type=str, choices=["PE", "SE"],help="sequencing type: PE or SE ? (default=PE)", default="PE")

	group2 = parser.add_argument_group('Advanced', 'optional arguments')
	group2.add_argument('--allele-frq', help='calculate allele frequency, if you just want calculate frequency only, please run it with "--run-split", and change the input file to "XXXX.reads.list" (default=False) ', action="store_true")
	group2.add_argument('--allele-frq-depth', help='run with "--allele-frq", the minimum reads depth in each snps which used for calculating allele frequency (default=20)', type=int, default=20)
	group2.add_argument('--bin-depth', help='the minimum reads depth in each windows (default=2000)', type=int, default=2000)
	group2.add_argument('--parallel-sort', help='Use GUN parallel for sorting', action="store_true")

	group3 = parser.add_argument_group('Sub-modules', 'run separately arguments')
	group3.add_argument('--run-split', help='run pipeline separately (default=False)', action="store_true")
	group3.add_argument('--step1', help='run step-1: Collecting mapped reads intersect with snps (by bedtools)', action="store_true")
	group3.add_argument('--step2', help='run step-2: Scan the snps and list the properly mapped reads in every snps', action="store_true")
	group3.add_argument('--step3', help='run step-3: Sort reads list file and filter out the reads which only contain 1 snp', action="store_true")
	group3.add_argument('--step4', help='run step-4: Identify the snps status in each reads', action="store_true")
	group3.add_argument('--step5', help='run step-5: Divide reads into non-recombination area and recombination area', action="store_true")
	group3.add_argument('--step6', help='run step-6: Sort the files and calculate the number of reads in different areas', action="store_true")
	group3.add_argument('--step8', help='run step-8: make windows, calculate the recombination reads number and non-recombination reads number in each bin', action="store_true")

	return parser

def quality_control(inputfile, snps, mapq, outputfile):
	print ("------ Quality control (by samtools) ------")
	time_start = time.time()

	if os.path.exists(outputfile + ".Q_" + str(mapq) + ".bam"):
		print (outputfile + ".Q_" + str(mapq) + ".bam existed, pass..")
		pass
	else:
		cmd1 = "samtools view -@ 24 -bq " + str(mapq) + " " + inputfile + " > " + outputfile + ".Q_" + str(mapq) + ".bam"
		print ("CMD = " + cmd1)
		os.system(cmd1)

	print ("------ END ------")
	print ("time: " + str (time.time()-time_start))

	return (outputfile + ".Q_" + str(mapq) + ".bam")

def step1(inputfile, snps, outputfile):
	print ("------ Step1: Collecting mapped reads intersect with snps (by bedtools) ------")
	time_start = time.time()

	if os.path.exists(outputfile + ".ref.snps.bed"):
		print (outputfile + ".ref.snps.bed existed, pass..")
		pass
	else:
		cmd1 = "awk '{print" + ' $1"\\t"$2"\\t"$2"\\t"$3"\\t"$4' + "}' " + snps + " |sort -nk 1 -nk 2 > " + outputfile + ".ref.snps.bed"
		print ("CMD = " + cmd1)
		os.system(cmd1)
	
	if os.path.exists(outputfile + ".filter.bam"):
		print (outputfile + ".filter.bam existed, pass..")
		pass
	else:
		cmd2 = "bedtools intersect -abam " + inputfile + " -b " + outputfile + ".ref.snps.bed" + " > " + outputfile + ".filter.bam"
		print ("CMD = " + cmd2)
		os.system(cmd2)

	print ("------ Step1 END ------")
	print ("time: " + str (time.time()-time_start))

	return (outputfile + ".filter.bam")

def step2(inputfile, snps, outputfile, tp, NM): # This step use the vcf file, should be fixed in the next version
	print ("------ Step2: Scan the snps and list the properly mapped reads in every snps ------")
	time_start = time.time()


	# build the index for bam file
	if os.path.exists(inputfile + ".bai"):
		print (inputfile + " index existed, pass..")
		pass
	else:
		cmd = "samtools index " + inputfile
		print ("CMD = " + cmd)
		os.system(cmd)


	fout = open(outputfile + ".reads.list", 'w')
	samfile = pysam.AlignmentFile(inputfile, "rb")
	
	if tp == "PE":
		with open(snps) as f:
			for each in f:
				each = each[:-1].split()
				#print (each)
				for pileupcolumn in samfile.pileup(each[0], int(each[1]), int(each[1]) + 1):
					#print (pileupcolumn)
					
					if pileupcolumn.pos != int(each[1]) - 1: 
						continue
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.get_tag("NM") > int(NM):
							continue
							print (pileupread.alignment)
						if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_proper_pair:
							#print (pileupread)
							fout.write (pileupread.alignment.query_name + "\t" + pileupread.alignment.query_sequence[pileupread.query_position] + "\t" + each[2] + "\t" + each[3] + "\t" + each[0] + "\t" + each[1] + "\t" + str(pileupread.alignment.query_qualities[pileupread.query_position]) + "\n")

	if tp == "SE":
		with open(snps) as f:
			for each in f:
				each = each[:-1].split("\t")
				#print (each)
				for pileupcolumn in samfile.pileup(each[0], int(each[1]), int(each[1]) + 1):
					if pileupcolumn.pos != int(each[1]) - 1: 
						continue
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.get_tag("NM") > int(NM):
							continue
						if not pileupread.is_del and not pileupread.is_refskip:
							fout.write (pileupread.alignment.query_name + "\t" + pileupread.alignment.query_sequence[pileupread.query_position] + "\t" + each[2] + "\t" + each[3] + "\t" + each[0] + "\t" + each[1] + "\t" + str(pileupread.alignment.query_qualities[pileupread.query_position]) + "\n")



	print ("------ Step2 END ------")
	print ("time: " + str (time.time()-time_start))

	return (outputfile + ".reads.list")
	
def calculate_frq(inputfile, snps, fai, outputfile, depth): # calculate allele frequence in each snps
	print ("------ Calculate frq: calculate allele frequence in each snps ------")
	time_start = time.time()

#	'''

	fout = open(outputfile + ".allele_frq", 'w')
	fout.write("chr\tpos\tallele_A\tallele_B\tfrq\n")

	with open(inputfile) as f:
		_, *ff, _e = f
		tmp = _[:-1].split("\t")
		chro_tmp = tmp[4]
		pos_tmp = tmp[5]

		if tmp[1] == tmp[2]:
			num_a = 1
			num_b = 0
		elif tmp[1] == tmp[3]:
			num_a = 0
			num_b = 1

		for i in ff:
			i = i[:-1].split("\t")
			chro = i[4]
			pos = i[5]
			if chro == chro_tmp and pos == pos_tmp:
				if i[1] == i[2]:
					num_a += 1
				elif i[1] == i[3]:
					num_b += 1
			else:
#				print ("--- " + )
				if num_a == 0 and num_b == 0:
					print ("chr " + chro_tmp + " pos " + pos_tmp + " is zero !!!")
				else:
					fout.write(chro_tmp + "\t" + pos_tmp + "\t" + str(num_a) + "\t" + str(num_b) + "\t" + str(num_a/(num_a + num_b)) + "\n")
#					print ("--- " + chro + " --- " + pos)
				chro_tmp = chro
				pos_tmp = pos
				num_a = 0
				num_b = 0
				if i[1] == i[2]:
					num_a += 1
				elif i[1] == i[3]:
					num_b += 1

		tmp = _e[:-1].split("\t")
		chro = tmp[4]
		pos = tmp[5]
		if chro == chro_tmp and pos == pos_tmp:
			if tmp[1] == tmp[2]:
				num_a += 1
			elif tmp[1] == tmp[3]:
				num_b += 1
			fout.write(chro_tmp + "\t" + pos_tmp + "\t" + str(num_a) + "\t" + str(num_b) + "\t" + str(num_a/(num_a + num_b)) + "\n")
		else:
			chro_tmp = chro
			pos_tmp = pos
			if tmp[1] == tmp[2]:
				num_a = 1
				num_b = 0
			elif tmp[1] == tmp[3]:
				num_a = 0
				num_b = 1
			fout.write(chro_tmp + "\t" + pos_tmp + "\t" + str(num_a) + "\t" + str(num_b) + "\t" + str(num_a/(num_a + num_b)) + "\n")

	fout.close()
	print ("done...make windows...")

#	'''

	### make windows
	cmd1 = "bedtools makewindows -g " + fai + " -w 1000000 |awk '{print" + ' $0"\\ta"NR' + "}'|sort -k1,1 -k2,2n > " + "genome.windows"
	print ("CMD = " + cmd1)
	os.system(cmd1)

	cmd2 = "awk '{if(($3+$4)>=" + str(depth) + "){print" + ' $1"\\t"$2"\\t"$2+1"\\t"$3"\\t"$4"\\t"$5' + "}}' " + outputfile + ".allele_frq" + " |bedtools map -a ./genome.windows -b - -c 4 |awk 'BEGIN" + '{print "chr\\tstart\\tend\\tid\\tnumber"}{if($5=="."){print $1"\\t"$2"\\t"$3"\\t"$4"\\t0"}' + "else{print}}' > allele_A.bin"
	print ("CMD = " + cmd2)
	os.system(cmd2)

	cmd3 = "awk '{if(($3+$4)>=" + str(depth) + "){print" + ' $1"\\t"$2"\\t"$2+1"\\t"$3"\\t"$4"\\t"$5' + "}}' " + outputfile + ".allele_frq" + " |bedtools map -a ./genome.windows -b - -c 5 |awk 'BEGIN" + '{print "chr\\tstart\\tend\\tid\\tnumber"}{if($5=="."){print $1"\\t"$2"\\t"$3"\\t"$4"\\t0"}' + "else{print}}' > allele_B.bin"
	print ("CMD = " + cmd3)
	os.system(cmd3)

	cmd4 = "paste allele_A.bin allele_B.bin |awk 'BEGIN{print" + ' "chr\\tstart\\tend\\tid\\tnumber"}{if(($5!="0")&&($1!="chr")){print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5/($5+$10)' + "}}' > " + outputfile + ".allele_frq.bin"
	print ("CMD = " + cmd4)
	os.system(cmd4)

	print ("------ END ------")
	print ("time: " + str (time.time()-time_start))
				


def step3(inputfile, snps, outputfile, qs, parallel):
	print ("------ Step3: Sort reads list file and filter out the reads which only contain 1 snp ------")
	time_start = time.time()

	if parallel:
		if int(qs) != 0:
			cmd = "awk '{if($7>=" + str(int(qs)) + "){print}}' " + inputfile + " | parallel --pipe --files sort -T $(pwd) -nk 1 -S512M | parallel -Xj1 sort -nk 1 -S1024M -m {} ';' rm {} |awk 'BEGIN{a=0}{if($1!=a){a=$1; tmp=$0; f=0}else{if(f==0){print tmp;print $0;f=1}else{print $0}}}' > " + outputfile + ".sorted.reads.list"
			print ("CMD = " + cmd)
			os.system(cmd)
		else:
			cmd = "cat " + inputfile + " | parallel --pipe --files sort -T $(pwd) -nk 1 -S512M | parallel -Xj1 sort -T $(pwd) -nk 1 -S1024M -m {} ';' rm {} |awk 'BEGIN{a=0}{if($1!=a){a=$1; tmp=$0; f=0}else{if(f==0){print tmp;print $0;f=1}else{print $0}}}' > " + outputfile + ".sorted.reads.list"
			print ("CMD = " + cmd)
			os.system(cmd)
	else:
		if int(qs) != 0:
			cmd = "awk '{if($7>=" + str(int(qs)) + "){print}}' " + inputfile + " | sort -T $(pwd) -nk 1 |awk 'BEGIN{a=0}{if($1!=a){a=$1; tmp=$0; f=0}else{if(f==0){print tmp;print $0;f=1}else{print $0}}}' > " + outputfile + ".sorted.reads.list"
			print ("CMD = " + cmd)
			os.system(cmd)
		else:
			cmd = "sort -T $(pwd) -nk 1 " + inputfile + " |awk 'BEGIN{a=0}{if($1!=a){a=$1; tmp=$0; f=0}else{if(f==0){print tmp;print $0;f=1}else{print $0}}}' > " + outputfile + ".sorted.reads.list"
			print ("CMD = " + cmd)
			os.system(cmd)

	print ("------ Step3 END ------")
	print ("time: " + str (time.time()-time_start))

	return (outputfile + ".sorted.reads.list")


def step4(inputfile, snps, outputfile):
	print ("------ Step4: Identify the snps status in each reads ------")
	time_start = time.time()

	fout = open(outputfile + ".sorted.reads.list.status", 'w')
	read = {}
	with open(inputfile) as f:
		for line in f:
			tmp = line.split("\t")[0]
			elements = line[:-1].split("\t")
			if tmp not in read:
				if len(list(read.keys())) > 0:
					status(fout, read)
				read = {}
				read[tmp] = [elements]
			else:
				read[tmp].append(elements)
				read[tmp] = sorted(read[tmp], key = lambda x:x[5]) 

	print ("------ Step4 END ------")
	print ("time: " + str (time.time()-time_start))

	return (outputfile + ".sorted.reads.list.status")

def status(fout, read):
	fout.write(list(read.keys())[0])
	for each in list(read.values())[0]:
		if each[1] == each[2]:
			status = "0"
		elif each[1] == each[3]:
			status = "1"
		else:#------------------------------------------ without status "2"
			#status = "2"
			continue
		fout.write("\t" + status + ":" + each[4] + ":" + each[5])
	fout.write("\n")

def step5(inputfile, snps, outputfile):
	print ("------ Step5: divide reads into non-recombination area and recombination area ------")
	time_start = time.time()

	f_no_rc = open(outputfile + ".no.recombination" , 'w')
	f_rc = open(outputfile + ".recombination" , 'w')

	with open(inputfile) as f:
		for i in f:
			i = i[:-1].split("\t")
			rc_num = 0
			rc_list = []
#			print ("---" + i[0])
			if len(i) >= 5:#---------------- add the judgment
#				print (i)
				#print (rc_list)
				start = i[1]
				length = 1
				for each in range(len(i) - 2):
					if i[each + 1].split(":")[0] == i[each + 2].split(":")[0]:
						length += 1
#					pass
					else:
						if each + 4 <= len(i):
							if length >= 2 and i[each + 2].split(":")[0] == i[each + 3].split(":")[0]:
								rc_list.append("\t".join(i[each + 1].split(":")[1:]) + "\t" + "\t".join(i[each + 2].split(":")[1:]) + "\t" + i[0] + "\n")
#							if start != i[each + 1]:
#								f_no_rc.write("\t".join(start.split(":")[1:]) + "\t" + "\t".join(i[each + 1].split(":")[1:]) + "\t" + i[0] + "\n")
								start = i[each + 2]
							else:
								rc_list.append("NULL")
								length = 1
								
						else:
							rc_list.append("NULL")
							break
						
#				if start != i[-1]:
				f_no_rc.write("\t".join(i[1].split(":")[1:]) + "\t" + "\t".join(i[-1].split(":")[1:]) + "\t" + i[0] + "\n")
				if len(rc_list) == 1 and rc_list[0] != "NULL":
#					print (rc_list)
					f_rc.write(rc_list[0])
				

	print ("------ Step5 END ------")
	print ("time: " + str (time.time()-time_start))

	return outputfile

def step6(inputfile, snps, outputfile, parallel):
	print ("------ Step6: sort the files and calculate the number of reads in different areas ------")
	time_start = time.time()

	if parallel:
		cmd1 = "cut -f1-4 " + inputfile + ".no.recombination" + " | parallel --pipe --files sort -k1,1 -k2,2n -S512M | parallel -Xj1 sort -k1,1 -k2,2n -S1024M -m {} ';' rm {} | awk '{if(NR==1){tmp=$0;a=0}if($0==tmp){a++;}" + 'else{print tmp"\\t"a; tmp=$0; a=1}}END{print tmp"\\t"a}' + "' > " + outputfile + ".no.recombination.number"
		cmd2 = "cut -f1-4 " + inputfile + ".recombination" + " | parallel --pipe --files sort -k1,1 -k2,2n -S512M | parallel -Xj1 sort -k1,1 -k2,2n -S1024M -m {} ';' rm {} | awk '{if(NR==1){tmp=$0;a=0}if($0==tmp){a++;}" + 'else{print tmp"\\t"a; tmp=$0; a=1}}END{print tmp"\\t"a}' + "' > " + outputfile + ".recombination.number"
	else:
		cmd1 = "cut -f1-4 " + inputfile + ".no.recombination" + " | sort -k1,1 -k2,2n | awk '{if(NR==1){tmp=$0;a=0}if($0==tmp){a++;}" + 'else{print tmp"\\t"a; tmp=$0; a=1}}END{print tmp"\\t"a}' + "' > " + outputfile + ".no.recombination.number"
		cmd2 = "cut -f1-4 " + inputfile + ".recombination" + " | sort -k1,1 -k2,2n | awk '{if(NR==1){tmp=$0;a=0}if($0==tmp){a++;}" + 'else{print tmp"\\t"a; tmp=$0; a=1}}END{print tmp"\\t"a}' + "' > " + outputfile + ".recombination.number"
	
	print ("CMD = " + cmd1)
	os.system(cmd1)

	print ("CMD = " + cmd2)
	os.system(cmd2)

	print ("------ Step6 END ------")
	print ("time: " + str (time.time()-time_start))

	return outputfile

def step7(inputfile, snps, outputfile):
	print ("------ Step7: for each recombination area, extract non-recombination reads number and recombination reads number, calculate rate ------")
	time_start = time.time()

	f_no_rc = open(outputfile + ".no.recombination.number" , 'r').read().split("\n")
	fout = open(outputfile + ".rate", 'w')
	with open(outputfile + ".recombination.number") as f:
		for i in f:
			rc = i[:-1].split("\t")
			rc_num = rc[-1]
			nrc_num = 0
			del_lines = []
			for nrc_line in f_no_rc:
				nrc = nrc_line.split("\t")

				### insert the judgment for distinguishing the chromosomes
				if nrc[0] == rc[0]: 


					if int(nrc[3]) < int(rc[1]): # no-rc.right < rc.left
						#print ("1---" + nrc_line)
						del_lines.append(nrc_line)
	#					f_no_rc.remove(nrc_line)
						#print ("1---rm---" + f_no_rc[0])
					elif int(nrc[1]) <= int(rc[1]) and int(nrc[3] >= rc[3]):
						#print ("2--" + nrc_line)
						nrc_num += int(nrc[-1])
					elif int(nrc[1]) > int(rc[3]): # no-rc.left > rc.right
						#print ("3--" + nrc_line)
						break
			for d in del_lines:
				f_no_rc.remove(d)
			#print (len(f_no_rc))
			fout.write (rc[0] + "\t" + rc[1] + "\t" + rc[2] + "\t" + rc[3] + "\t" + rc_num + "\t" + str(nrc_num) + "\t" + str(int(rc_num)/(int(nrc_num) + int(rc_num))) + "\n")
	#print (f_no_rc)

	print ("------ Step7 END ------")
	print ("time: " + str (time.time()-time_start))

	return outputfile

def step8(inputfile, fai, outputfile, win, step, bin_depth, tp, insertion):
	print ("------ Step8: make windows, calculate the recombination reads number and non-recombination reads number in each bin ------")
	time_start = time.time()

	cmd1 = "bedtools makewindows -g " + fai + " -w " + str(win) + " -s " + str(step) + " |awk '{print" + ' $0"\\ta"NR' + "}'|sort -k1,1 -k2,2n > " + "genome.windows"
	print ("CMD = " + cmd1)
	os.system(cmd1)
	
	cmd2 = "awk '{if" + '($2<$4){print $1"\\t"$2"\\t"$4"\\t"$5}' + \
		"}' " + outputfile + ".no.recombination.number | sort -k1,1 -k2,2n | " + "bedtools map -a ./genome.windows -b - -c 4 |" + \
		"awk 'BEGIN{print" + ' "chr\\tstart\\tend\\tid\\tnumber"}{if($5=="."){print $1"\\t"$2"\\t"$3"\\t"$4"\\t0"}' + "else{print}}' >" + outputfile +  ".no.recombination.bin"
	print ("CMD = " + cmd2)
	os.system(cmd2)
	
	cmd3 = "awk '{if" + '($2<$4){print $1"\\t"$2"\\t"$4"\\t"$5}' + \
		"}' " + outputfile + ".recombination.number | sort -k1,1 -k2,2n | " + "bedtools map -a ./genome.windows -b - -c 4 |" + \
		"awk 'BEGIN{print" + ' "chr\\tstart\\tend\\tid\\tnumber"}{if($5=="."){print $1"\\t"$2"\\t"$3"\\t"$4"\\t0"}' + "else{print}}' >" + outputfile +  ".recombination.bin"
	print ("CMD = " + cmd3)
	os.system(cmd3)

	if tp == "PE":
		cmd4 = "paste " + outputfile +  ".recombination.bin " + outputfile +  ".no.recombination.bin" + \
			" | awk 'BEGIN{print" + ' "chr\\tstart\\tend\\tid\\tnumber"}{if(($5!="0")&&($1!="chr")&&($10>=' + str(bin_depth) + ')){print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5/($5+$10)/' + str(insertion) + r'*' + str(win) + "}}' >" + outputfile + ".bin.rate"
	else:
		cmd4 = "paste " + outputfile +  ".recombination.bin " + outputfile +  ".no.recombination.bin" + \
			" | awk 'BEGIN{print" + ' "chr\\tstart\\tend\\tid\\tnumber"}{if(($5!="0")&&($1!="chr")&&($10>=' + str(bin_depth) + ')){print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5/($5+$10)}' + "}' >" + outputfile + ".bin.rate"
	print ("CMD = " + cmd4)
	os.system(cmd4)
	
	print ("------ Step8 END ------")
	print ("time: " + str (time.time()-time_start))

	return outputfile


def main():
	parser = get_parser()
	args = vars(parser.parse_args())

	if args['input'] is not None:
		print(version())
	
	inputfile = args['input']

	if args['run_split']:
		print ("------- Run the pipeline separately -------")
		if args['mapq'] == 0:
			pass
		else:   
			inputfile = quality_control(inputfile, args['snps'], args['mapq'], args['output'])
		if args['step1']:
			inputfile = step1(inputfile, args['snps'], args['output'])
		if args['step2']:
			inputfile = step2(inputfile, args['snps'], args['output'], args['type'], args['NM'])
		if args['allele_frq']:
			calculate_frq(inputfile, args['snps'], args['fai'], args['output'], args['allele_frq_depth'])
		if args['step3']:	
			inputfile = step3(inputfile, args['snps'], args['output'], args['qs'], args['parallel_sort'])
		if args['step4']:
			inputfile = step4(inputfile, args['snps'], args['output'])
		if args['step5']:
			inputfile = step5(inputfile, args['snps'], args['output'])
		if args['step6']:
			inputfile = step6(inputfile, args['snps'], args['output'], args['parallel_sort'])
		#	inputfile = step7(inputfile, args['snps'], args['output'])
		if args['step8']:
			inputfile = step8(inputfile, args['fai'], args['output'], args['window'], args['step'], args['bin_depth'], args['tp'], args['insertion'])
		print ("------- Done!!! -------")
	
	else:
		if args['mapq'] == 0:
			pass
		else:
			inputfile = quality_control(inputfile, args['snps'], args['mapq'], args['output'])

		inputfile = step1(inputfile, args['snps'], args['output'])
		inputfile = step2(inputfile, args['snps'], args['output'], args['type'], args['NM'])

		if args['allele_frq']:
			calculate_frq(inputfile, args['snps'], args['fai'], args['output'], args['allele_frq_depth'])
	
		inputfile = step3(inputfile, args['snps'], args['output'], args['qs'], args['parallel_sort'])
		inputfile = step4(inputfile, args['snps'], args['output'])
		inputfile = step5(inputfile, args['snps'], args['output'])
		inputfile = step6(inputfile, args['snps'], args['output'], args['parallel_sort'])
		#	inputfile = step7(inputfile, args['snps'], args['output'])
		inputfile = step8(inputfile, args['fai'], args['output'], args['window'], args['step'], args['bin_depth'], args['tp'], args['insertion'])


if __name__ == "__main__":
	main()
