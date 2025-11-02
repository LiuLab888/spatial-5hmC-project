#######################
General description                    
#######################

All details in the spatial 5hmC analysis process are outlined in this file, which includes data preprocessing, read alignment, hydroxymethylation extraction, spatial clustering of hydroxymethylation, DhMR identification, and GO enrichment analysis.


#######################
Spatial hydroxymethylation analysis pipeline
#######################

All codes in Spatial_methyl_pipeline/spatial_template_cmd_scaled.ibp.sh contains steps of data preprocessing, read alignment, and hydroxymethylation extraction, while the files in /bins and /python_script are invoked scripts for supporting the normal operation of the core script spatial_template_cmd_scaled.ibp.sh. It should be emphasized that the user needs to substitute the sample name, fastq file path, the number of pixels, and the spatial location list with those you use in the corresponding brace, as shown below.

sed -e "s/sample_name/${name}/g" spatial_template_cmd_scaled.ibp.sh > ${name}_spatial_template_cmd_scaled.ibp.sh
sed -i "s/\/dir/\/${dir}/g" ${name}_spatial_template_cmd_scaled.ibp.sh
sed -i "s/all_cells/${all_cells}/g" ${name}_spatial_template_cmd_scaled.ibp.sh
sed -i "s/each_part/${each_part}/g" ${name}_spatial_template_cmd_scaled.ibp.sh
sed -i "s/combine_barcode.round2round1_index1_index2.v3_big_hydroxymethylation.txt/${dir}_ref_barcode.txt/g" ${name}_spatial_template_cmd_scaled.ibp.sh

Eventually, you should run the following command on a Linux system.

sh ${name}_spatial_template_cmd_scaled.ibp.sh

After these steps are completed, you will have all the preliminary results for downstream analysis, including the raw read matrix in A_prework_outputdir_${sample_name}/fastq, the coverage, depth, and hydroxymethylation matrices in A_prework_outputdir_${sample_name}/bismark/split/, and the hydroxymethylation bed files of all pixels in A_prework_outputdir_${sample_name}/bismark/split/part*.  

If you experience any environment configuration error, please contact us at liujiang@ibp.ac.cn


#######################
Spatial hydroxymethylation clustering
#######################

Owing to the novel spatial clustering method we designed, it is not well-packaged, which means that many flexible parameters need to be tuned by the user for the best clustering performance. Herein, we give a clear clustering process and the tools used, convenient for everyone to build their spatial hydroxymethylation clustering code. We also exhibit all clustering scripts with adjusted parameters shown in our study in Spatial_methyl_clustering/spatial_hydroxymethylation_cluster.R and the folder /data contains all related data. All steps of spatial clustering are as follows.

1.	Obtain the hydroxymethylation level matrix in A_prework_outputdir_${sample_name}/bismark/split/
2.	Build the Euclidean distance matrix for all pixels (e.g., for a slice with 96 pixels*96 pixels, a matrix with 9216 rows and 9216 columns)
3.	Identify the critical anchors using the FindIntegrationAnchors function in Seurat to construct the distance matrix after dimensionality reduction
4.	Weight the dimensionally reduced spatial distance matrix and the hydroxymethylation level matrix with a series of weight coefficients (Suggest using equal weights) and merge the two by column.
5.	Perform the K-means combined with get_clust_tendency and clusGap function in the R package factoextra v1.0.7 to evaluate the clustering tendency and the gap statistic for identifying the optimal cluster number
6.	Carry out the Uniform Manifold Approximation and Projection (UMAP) analysis for data visualization, after determining all principal components in PCA analysis (if fewer columns, the PCA analysis can be skipped). 
7.	Project all clusters to the spatial location and visualize the spatial cluster pattern using the ggplot2 function in R package.

If you need any detailed discussion for the above clustering method, please send an email to liujiang@ibp.ac.cn


#######################
Co-clustering analysis
#######################

All codes in /Co_clustering_analysis/Co_clustering_analysis.R show the analysis of unsupervised co-clustering of our SmC-seq data with previously published single-cell hydroxymethylation data from mouse embryos. 

The user needs to generate the hydroxymethylation level and CpG coverage number of all genebody regions in advance for all pixels or cells, and then follow the step of data input in the script to generate the two matrices for data quality control and subsequent clustering analysis. To display the feasibility of our code, we provide the two matrices we have generated in /data (Due to the large storage space requirement, the two matrices are headed with 1000 rows).


#######################
DhMR analysis
#######################

In DhMR analysis, the user needs to merge the bed files for all pixels with the same cluster using DhMR_analysis/01_merge_bismark_bed_files.sh, and then transform the format of the merged bed files into the standard input by DhMR_analysis/02_pre_DhMR_analysis.sh, followed by the use of callDhMR function in the package DSS for automatically identifying DhMRs or utilizing the fisher.test function in the package stats for custom-sized DhMRs in DhMR_analysis/03_DhMR_analysis.R.


#######################
GO enrichment analysis
#######################

We established a routine analysis process in GO_analysis/GO_analysis.R using the clusterProfiler package and also provided some peak data in /data for convenient validation of the R script's usability.


#######################
Spatial RNA analysis
#######################
Due to the availability of comprehensive spatial RNA analysis methods currently, we perform the ST pipeline, which is widely used to perform data processing at spatial transcriptomics. Our similar processes are shown in spatial_RNA_analysis/, and more details about the developer can be found at https://github.com/jfnavarro/st_pipeline. 

