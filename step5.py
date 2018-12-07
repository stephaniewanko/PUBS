import pandas as pd
import argparse

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="find change in fitness values between treated and untreated")
	required = parser.add_argument_group('required')
	required.add_argument('-t', '--csv_treated',  nargs='*',required=True,
		help='fitness values of treated condition')
	required.add_argument('-u', '--csv_untreated',  nargs='*',required=True,
		help='fitness values of untreated condition')
	required.add_argument('-n', '--name',  nargs='*',required=True,
		help='name of compound to label output csv file')
	args = parser.parse_args()
	ct = pd.read_csv(args.csv_treated[0],index_col='Unnamed: 0')
	cu = pd.read_csv(args.csv_untreated[0],index_col='Unnamed: 0')
	mergedF = pd.merge(ct, cu, on='variant')
	mergedF['changeInFitness'] = mergedF['fitness_x'] - mergedF['fitness_y']
	mergedF.to_csv(args.name[0] + "barcode_fitness_change.csv")

