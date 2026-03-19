#!/usr/bin/env python

import argparse
import subprocess
import os.path
import gzip 

def info_db(db_info_file):
	name=db_info_file.split()[0]
	path =os.path.realpath(db_info_file.split()[1])

	if db_info_file.split()[1]=='nodes.dmp':
		path=os.path.dirname(path)

	# get db size 
	process = subprocess.run(['du', '-sh',path], stdout=subprocess.PIPE)
	size = process.stdout.split()[0].decode('utf-8')

	if (name == "Host_genome"):
		size = f"{size} ({get_genome_seq_count(path)} seq)"
	
	# get date of last modifaction 
	process = subprocess.run(['stat', '-c %y', path], stdout=subprocess.PIPE)
	modif_date = process.stdout.split()[0].decode('utf-8')

	return name,size,modif_date,path


def get_genome_seq_count(genome_path):

	seq_count = 0

	proper_open = gzip.open if genome_path.endswith('.gz') else open
	with proper_open(genome_path,"rt") as fh:
		for l in fh:
			if l.startswith('>'):
				seq_count += 1 
	return seq_count 

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--all_db', required = True, help = \
	'text file with details of all databank used ')
	args = parser.parse_args()
	
	print("DB (folder or file)\tSize\tLast modification\tPath")
	
	with open(args.all_db,"r") as list_db:
		lines = list_db.readlines()
		for line in lines:
			print(line)
			list_info=info_db(line.strip())
			print(*list_info,sep="\t")



if __name__ == "__main__":
    main()
