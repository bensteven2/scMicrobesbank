# scMicrobesbank
#1 Microbial_Data_Extraction.sh
#This script is designed to extract microbiome data from sequencing data. It performs several steps of analysis and utilizes various tools and software packages. 
#Usage:
#sh scripts_nature.sh sample_name mem fqPath analysis_folder
#Parameters:
#sample_name: The name of the sample.
#mem: Maximum memory usage.
#fqPath: The path of the sample_fq.
#analysis_folder: The path of the analysis.
#2.1 Data_Quality_Control_And_Cell_Annotation_Example _GSE166635.ipynb
#This script is designed to perform data quality control and batch effect correction on single-cell RNA sequencing (scRNA-seq) data. The script utilizes the scanpy #library (v1.8.2) and implements various techniques for filtering low-quality cells, removing rarely expressed genes, normalizing the data to UPM (UMI per million), #integrating across samples, and removing batch effects using Harmony and BBKNN.Cell clustering is performed using the Leiden algorithm based on the graph.Cell type annotation is processed based on the methods described in the corresponding paper. If available, metadata provided by the original authors are used for cell type annotations.
#2.2 Processing_Microbiome_File_Example_GS
#Visualization_of_Downstream_Analysis_of_Microbiome_Data_Example_GSE166635.ipynb
#Nebula.R
