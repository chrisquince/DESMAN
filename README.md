# DESMAN
De novo Extraction of Strains from MetAgeNomes

First step is to identify variant positions. The '-p' flag uses one dimenisonal optimisition to find individual base frequencies.

	python ./desman/Variant_Filter.py data/Mock_15_var_freq.csv -o COG0015_out -p

Generates output files: 

