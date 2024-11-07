# Isoform_quantification

To try the initial algorithm, please run the file 'AT_CODE/main_EM_VI_GD_CLEAN.py'. This is the main file and it will call EM_VI_GD_CLEAN.py and DirichletOptimizer_CLEAN.py. <br />

**DATA:** Please find some trial SIRV data in this folder: 'SIRV_data' <br />
<br />
**ENVIRONMENT:** Please install the environment from the following .yml or .txt file: 'AT_CODE/NanoCount_5.yml', 'AT_CODE/NanoCount_5.txt'<br />
<br />
**Command:** To run with trial SIRV data, enter this python command: <br />
conda activate NanoCount_5 <br />
python main_EM_VI_GD_CLEAN.py --input_folder "data folder name"  > "output file name" <br />
<br />
eg.python main_EM_VI_GD_CLEAN.py --input_folder /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln_pklfiles/ > output.txt <br />
<br />
the command will automatically create a folder called 'result' and store the isoform abundances. 
