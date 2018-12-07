import pickle
import csv
pickleFile = "barcode_to_variant_dict_cutoff_2.pkl"
pickleDict = pickle.load( open(pickleFile, "rb" ) )
with open('barcode_dict.csv', 'w') as csv_file:
	writer = csv.writer(csv_file)
	for key, value in pickleDict.items():
		pos, variant = value
		writer.writerow([key, pos, variant])
