from Bio import SeqIO
from collections import defaultdict
import csv
from multiprocessing import Pool
import argparse
import pickle

# input: sample fastqfile, number of bases to trim from back of barcode
# output: valid barcode dictionary, bad barcode dictionary
# TODO: recover mismatched with no matching barcode
def fastqToDict(fastqFile, barcodeLength, barcodeL, trim=True):
	fastqDict = defaultdict(int)
	badCodesDict = defaultdict(int)
	nCount = 0
	missingBarcodeCount = 0
	correctCount = 0
	for record in SeqIO.parse(fastqFile, "fastq"):
		currentCode = str(record.seq)
		if trim:
			trimLength = len(currentCode) - barcodeLength
			trimmedCode = currentCode[0:len(currentCode)-trimLength]
		else:
			trimmedCode = currentCode
		if ("N" in trimmedCode): 
			badCodesDict[str(trimmedCode)] += 1
			nCount += 1
		elif (trimmedCode in barcodeL):
		#else:
			fastqDict[str(trimmedCode)] += 1
			correctCount += 1
		else: # barcode not in barcode dictionary
			badCodesDict[str(trimmedCode)] += 1
			missingBarcodeCount += 1
	return (fastqDict, badCodesDict, nCount, missingBarcodeCount, correctCount)
	#return (nCount, missingBarcodeCount, correctCount)

# input: sample fastqfile, number of bases to trim from back of barcode
# output: valid barcode dictionary, bad barcode dictionary
# TODO: recover mismatched with no matching barcode
def fastqToDictRaw(fastqFile):
	fastqDict = defaultdict(int)
	for record in SeqIO.parse(fastqFile, "fastq"):
		currentCode = str(record.seq)[0:20]
		fastqDict[currentCode] += 1
	return (fastqDict)


# input: csv pickle file
# output: list of all valid barcodes in barcode dictionary
def csvToList(pickleCSV):
	barcodes = []
	with open(pickleCSV) as csv_file:
		csv_reader=csv.reader(csv_file)
		for row in csv_reader:
			barcodes.append(row[0])
	return barcodes

# input: two strings 
# output: hamming distance
def hamming(instrings):
	(s1, s2) = instrings
	dist=0
	for i in range(len(s1)):
		if s1[i] != s2[i]:
			dist += 1
		if dist > 1:
			return 2 # STOP if mismatch by more than 1
	return dist

def filterFun(input):
	(hamming, index) = input
	return hamming <= 1

# input: dictionary of barcodes and counts, minimum count to recover
# output: new, shorter badCodesDict
def filterBadCodes(badCodesDict, minCount):
	return {k: v for k, v in badCodesDict.items() if v >= minCount}

# input: dictionary of bad barcodes, list of valid barcodes
# output: list of barcodes to update and count to add
def recoverBadCodes(badCodesDict, barcodeL):
	newAddL = []
	pool = Pool(processes=4)
	recoveredCount = 0
	unrecoveredCount = 0
	for item in list(badCodesDict.items()):
		(badCode, count) = item
		badL = [badCode] * len(barcodeL)
		zippedL = list(zip(badL, barcodeL))
		mappedL = list(pool.map(hamming, zippedL))
		indicesL = list(range(len(mappedL)))
		totalMap = list(zip(mappedL, indicesL))
		filteredMap = list(filter(filterFun, totalMap))
		if len(filteredMap) == 1:
			recoveredCount += count
			(hammingNo, index) = filteredMap[0]
			goodCode = barcodeL[index]
			newAdd = (goodCode, count)
			newAddL.append(newAdd)	
		else:
			unrecoveredCount += count
	pool.close()	
	return (newAddL)

# input: countDict of barcodes/counts, list of updates to make to dictionary
# output: updated countDict
def updateCountDict(countDict,newUpdatesL):
	for update in newUpdatesL:
		(goodCode, count) = update
		countDict[goodCode] += count
	return countDict

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

# input: fastq file to count, barcode csv path
# output: nothing, write updated counts to new csv
def recoverNs(fastqFile, pickleDict, newCSVname):
	#barcodeL = csvToList(pickleCSV)
	barcodeL = list(pickleDict.keys())
	barcodeLength = len(barcodeL[0])
	barcodeS = set(barcodeL)
	(fastqDict, badCodesDict, nCount, missingBarcodeCount, correctCount) = fastqToDict(fastqFile, barcodeLength,barcodeS)
	print(nCount, missingBarcodeCount, correctCount)
	badCodesDict = filterBadCodes(badCodesDict, 1)
	newAddL = recoverBadCodes(badCodesDict, barcodeL)
	updatedDict = updateCountDict(fastqDict, newAddL)
	dictToCSV(updatedDict,newCSVname)
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="script to recover counts from barcodes that aren't in the dictionary")
	required = parser.add_argument_group('required')
	required.add_argument('-f', '--fastq',  nargs='*',required=True,
		help='raw fastq file with counts to recover')
	#required.add_argument('-b', '--barcode_pickle', required=True,
	#	help="pickle containing a dictionary with barcode as keys and library variants as values")
        required.add_argument('-b', '--barcode_csv', required=True,
                        help="csv containing a dictionary with barcode as keys and library variants as values")
        args = parser.parse_args()
        fastq = args.fastq[0]
	#with open(args.barcode_pickle, 'rb') as f:
        #        barcode_variant_dict = pickle.load(f)
        #barcodeL = list(barcode_variant_dict.keys())
        barcodeL = csvToList(args.barcode_csv)
        barcodeLength = len(barcodeL)
	barcodeS = set(barcodeL)
        print("Loaded dictionary")
	(fastqDict, badCodesDict, nCount, missingBarcodeCount, correctCount) = fastqToDict(fastq, barcodeLength,barcodeS)
	totalCount = nCount + missingBarcodeCount + correctCount
	missingCount = missingBarcodeCount + nCount
	print(str((missingCount*100.00)/totalCount) + "% of reads missing")
	badCodesDict = filterBadCodes(badCodesDict, 1)
	newAddL = recoverBadCodes(badCodesDict, barcodeL)
	updatedDict = updateCountDict(fastqDict, newAddL)
	newCSVname = ((args.fastq).split(".")[0]) + "_recovered.csv"
	dictToCSV(updatedDict,newCSVname)
        
