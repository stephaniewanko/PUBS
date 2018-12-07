from Bio import SeqIO
from collections import defaultdict
import csv
from multiprocessing import Pool
import argparse
import os

# input: sample fastqfile, number of bases to trim from back of barcode
# output: valid barcode dictionary, bad barcode dictionary
# TODO: recover mismatched with no matching barcode
def fastqToDictRaw(fastqFile):
	fastqDict = defaultdict(int)
	for record in SeqIO.parse(fastqFile, "fastq"):
		currentCode = str(record.seq)[0:20]
		fastqDict[currentCode] += 1
	return (fastqDict)

# input: count dictionary, name of new csv to write to
# output: None, writes to csv file
def dictToCSV(mydict,csvName):
	with open(csvName, 'w') as csv_file:
		writer = csv.writer(csv_file)
		header = ['barcode', 'count']
		writer.writerow(header)
		for key, value in mydict.items():
			writer.writerow([key, value])
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="script to write fastq to count csv")
	required = parser.add_argument_group('required')
	required.add_argument('-f', '--fastq',  nargs='*',required=True,
		help='raw fastq file with counts to recover')
	args = parser.parse_args()
        fastq = args.fastq[0]
        os.chdir("/data/lbny/recoveredReads/")
	fastqDict = fastqToDictRaw(fastq)
	newCSVname = ((fastq).split(".")[0]) + "_counts.csv"
	dictToCSV(fastqDict,newCSVname)
