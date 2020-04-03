ClonalTREE2

Requirements:
- Interpreters: Python 3, Perl
- Python libraries: numpy, graphviz
- Trimmomatic 0.33

Installation: 
git clone https://github.com/COL-IU/ClonalTREE2.git

Running the program: 

Part 1 - Aligning the FASTQ reads, variant calling and compilation of information in ClonalTREE2 recognizable formats*:
1) Make sure the FASTQ files are gzipped and named in the following format: projectname_timetag_1.fq.gz, projectname_timetag_2.fq.gz
2) Make changes to the file path variables and the timetags array appropriately in pipeline.sh
3) Run: pipeline.sh projectname

Part 2 - Running ClonalTREE2
1) Run: python3 ClonalTREE.py <VAF file> <rd file> <out prefix>\
VAF file:	[String] Input file containing the variant allele frequencies matrix (F).\
RD file:	[String] Input file containing the read depth matrix (R).\
out prefix:	[String] File path to prefix all output files.
  
File formats:
1) VAF file - variant allele frequencies matrix (F): First row should be space separated names of mutations that are column headers of the F matrix, which also will appear as node names in the inferred clonal tree. The subsequent rows should contain the allele frequencing of each mutation in the same order as the column headers, and the rows ordered by time points. 

Example VAF file:

1	2	3	4	5	6	7	8	9	10	11	12	13	14\
0.53850	0.51360	0.45410	0.47870	0.62280	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000\
0.43750	0.30560	0.48840	0.26190	0.34430	0.21880	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000\
0.41940	0.35850	0.64000	0.26830	0.28950	0.50000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000	0.00000\
0.58260	0.37500	0.63950	0.26390	0.31200	0.46430	0.07080	0.08860	0.20530	0.12290	0.20000	0.00000	0.00000	0.00000\
0.48210	0.49020	0.45710	0.75000	0.44830	0.60610	0.13510	0.17860	0.51850	0.39390	0.38600	0.10530	0.15380	0.61190\
0.55800	0.33720	0.63830	0.40540	0.36690	0.52520	0.00000	0.00000	0.28170	0.31360	0.33020	0.05190	0.25700	0.29960

2) RD file - read depth matrix (R): First row is same as the VAF file. Subsequent rows should contain the read depths at corresponding loci. 

Example rd file:

1	2	3	4	5	6	7	8	9	10	11	12	13	14\
182	220	207	188	289	0	0	0	0	0	0	0	0	0\
32	36	43	42	61	32	0	0	0	0	0	0	0	0\
31	53	50	41	38	44	0	0	0	0	0	0	0	0\
230	32	147	72	250	196	234	158	191	179	80	0	0	0\
56	51	35	36	58	33	77	56	54	34	57	57	52	67\
181	261	235	259	278	278	0	0	213	169	315	231	250	267

Note: These files are generated as a result of pipeline.sh

Output files:
1) <out prefix>.tree - Lists the tree with each node (clone) and corresponding ancestor. 
2) <out prefix>.dot - Dot source file to be used to visualize using graphviz.
3) <out prefix>.info - Additional useful information about the prediction. 
