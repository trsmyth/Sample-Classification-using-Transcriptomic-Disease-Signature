# Sample-Classification-using-Transcriptomic-Disease-Signature

This repository contains code for a personal side project seeking to classify human tissue samples according to transcriptomic signature based on publicly available disease labels (DiSignAtlas) using Medical Subject Heading (MeSH) designations. Overall, nothing was expected from this analysis, and it was primarily conducted to satisfy personal curiosity and to learn more about deep neural net (DNN) design and construction using PyTorch. 

In short, much more data is likely required to produce a DNN that is capable of reliably delineating disease signatures between closely related MeSH terms, if at all possible. This at least partially stems from overlap between related diseases and the MeSH terms which describe them, frequently leading to samples appearing in multiple term datasets. For example, [Asthma], [Lung diseases – obstructive], and [Hypersensitivity – immediate] are three MeSH terms which describe closely related situations and may never be differentiated, even if more samples are sorted into any of the specific categories. Collapsing highly related terms, potentially through an assessment of the fraction of overlapping samples, may be required. Further, as of now, a Random Forest approach (right) provides similar classification accuracy while requiring a fraction of the computational resources compared to a DNN (left). 

<img width="986" height="504" alt="image" src="https://github.com/user-attachments/assets/bc9e4fd4-5407-48cf-8351-dbc313aa796f" />

If this project was to be seriously developed, a regularization and feature selection approach such as Elastic Net would be more valuable. Determining key genes or groups of genes which correspond to MeSH terms, as well as the overlap between related terms, would provide more details regarding which MeSH terms should be collapsed as well as key information about important genes which may drive distinct or overlapping diseases. 

Brief description of methodology:
1.	Gene expression omnibus (GEO) sample designations corresponding to specific diseases and appropriate controls were collected from the disease signature atlas (DSA).
2.	Count data corresponding to these samples was isolated from the ARCHS4 RNAseq database, release v2.5.
3.	Diseases were assigned to the closest related medical subject heading (MeSH) designation, and MeSH tree branch points corresponding to those designations were collected.
4.	Samples were assigned a tissue designation (i.e. blood, airway epithelial cell, etc.) based on key word filtering as appears in GEO metadata.
5.	Samples were grouped according to their tissue type, MeSH term, and MeSH tree position to form disease signature sample groups.
6.	Count data was preprocessed, including gene and sample selection as well as batch correction.
7.	A deep neural net (DNN) was trained using the processed count data and metadata through a PyTorch workflow.
8.	A random forest model was grown using the ranger R library to compare the accuracy of a DNN (computationally expensive) vs an optimized random forest approach.

A detailed overview of the methodologies employed in this analysis and a description of the code presented in each folder:
1. Disease Signature Atlas:
  * Metadata for studies associated with various MeSH designations were downloaded from http://www.inbirg.com/disignatlas/browse -> Browse by disease tree -> All Diseases -> MeSH. This provides associated study IDs, disease names, tissue types, data type, and organism.
  * ‘Disease Information of All DiSignAtlas Datasets’ was downloaded from http://www.inbirg.com/disignatlas/download. This provides the study ID and the number of samples which are controls versus diseased in each dataset.
  * Rselenium and Rvest were used to scrape Control/Disease sample designations from http://www.inbirg.com/disignatlas/detail/{studyID}. This is necessary as sample gene expression omnibus (GEO) names are not always in sequential order related to X control/Y diseased samples.
  * Samples were given control or disease designations based on the number of listed control/diseased samples in ‘Disease Information of All DiSignAtlas Datasets’ and the order of appearance as found from scraping.
2.	ARCHS4:
  * The ARCHS4 database (version 2.5, accessed 11/23/2024) was downloaded from https://maayanlab.cloud/archs4/download.html as an h5 file.
  * Count data for samples found in both the DSA and ARCHS4 dataset were isolated and saved as a gene (row) by sample (column) matrix.
    * This analysis was conducted on my personal desktop which was unable to process the large volume of data at once. Instead, samples were extracted in 2500 sample chunks and merged later.
3.	MeSH Assignment:
  * Disease names associated with individual DSA studies were output as a .csv file and were aligned with the closest matching MeSH term through the NCBI MeSH database at https://www.ncbi.nlm.nih.gov/mesh and/or general term similarities. For example, Asthma, Allergric Asthma, and Eosinophilic Asthma were all assigned ‘Asthma’.
    * Note: While this was chosen in an attempt to find a pan-disease signature (such as Asthma in this example), it is very possible that no such signature exists or is very weak for a given disease when considering specific subsets. For example, Allergic Asthma may actually be better lumped with hypersensitivity while dampening a ‘base’ Asthma signature. An initial similarity scoring could be considered here to provide a clearer rationale for term matching.
  * Every associated MeSH tree position for the designated MeSH terms were isolated. This generally led to 1-4, but as many as 8, tree positions.
4.	Sample Tissue Designation:
  * The source tissue of sample type was assigned using a key word filtering approach.
    * I attempted to find the most distilled terms which applied as broadly as possible, but slight variations in naming conventions between groups, spelling errors, and large variations in naming details frequently forced the use of very specific or unique key words for a small number of studies or samples.
    * My approaches to sample designations have gotten much more refined in newer projects since conducting this phase. If I were to start over, a better approach would be to merge all metadata sections into a single metadata string followed by standardizing letter case and removal of punctuation. This, plus the addition of word boundaries, would have cut down the effort for these assignments drastically.
5.	Metadata Alignment:
  * Samples related to pulmonary diseases and control samples associated with tissues appearing in those datasets were isolated.
  * For diseased samples, all previously isolated tree positions associated with assigned MeSH terms were pruned to the third branch.
    * Control samples were designated ‘Control – [Tissue]’ and shared between diseased samples. In other words, a healthy control blood sample from every study which used blood were combined, regardless of whether that sample served as a control for disease x or disease y.
  * MeSH tree terms and their corresponding branch positions were downloaded from https://nlmpubs.nlm.nih.gov/projects/mesh/2024/meshtrees/ and MeSH terms associated with each new pruned branch position were assigned. For example, a sample which was described as having come from an asthmatic individual was assigned:
    * Asthma [C08.127.108] -> Asthma [C08.127.108]
    * Asthma [C08.381.495.108] -> Lung Diseases, Obstructive [C08.381.495]
    * Asthma [C08.674.095] -> Asthma [C08.674.095]
    * Asthma [C20.543.480.680.095] -> Hypersensitivity, Immediate [C20.543.480]
  * Samples were grouped according to their new MeSH designations and groups with less than 10 samples were removed from analysis.
    * In the example above, every Asthma sample would be included in three groups: [Asthma], [Lung Diseases, Obstructive], and [Hypersensitivity, Immediate].
6.	Count Data Preprocessing:
  * Metadata was merged with count data, duplicating samples which appear in multiple MeSH designation groups.
  * Data was converted to log count per million (LCPM) and genes were selected for inclusion using the edgeR function filterByExpr() followed by selection based on overall and within group variances.
    * FilterByExpr() provides an algorithmic approach to gene inclusion/exclusion rather than a basic cutoff like a mean count above a set threshold. See https://rdrr.io/bioc/edgeR/man/filterByExpr.html for further details.
    * Genes with an average within group variance of less than 5 and a ratio of overall variance to group variance above 5 (5x more overall variance than average within group variance) were retained.
  * Batch correction was conducted using the limma function removeBatchEffect() for both control and diseased samples.
    * See https://rdrr.io/bioc/limma/man/removeBatchEffect.html for further details.
    * Experimental series (GEO series ID) was treated as the batch term for control samples.
    * Experimental series AND tissue type were treated as batch and batch2, respectively, for diseased samples.
  This approach is absolutely off label. The desire is to create a ‘tissue-less’ gene signature for a given disease, but MUCH more research would have to be conducted into the statistical underlying of this approach to be comfortable using it in anything beyond a personal learning project.
  * Sample outliers were removed from the dataset using a dendrogram clustering approach.
  * Sample correlations were assessed for each group and dendrograms were assembled with a cut height of 0.25 (corresponding to Pearson correlation of 0.75).
  * Samples were split into a training and testing dataset using a 70/30 split for each group and merged with its corresponding metadata.
    * Note: This means each group was split 70/30 rather than a 70/30 split for the total dataset. This was done to ensure extreme minority classes were not unevenly split or entirely placed into either dataset.
7.	ARCHS4 DNN:
  * A deep neural net (DNN) was constructed and trained over 500 epochs.
  * Count data was converted to half precision (floating point 16) before training to save memory.
  * Cross entropy loss with balanced class weights and mean reduction was employed as the loss criteria.
    * Overall, with the current form of this data, the design of the DNN (number of layers and number of neurons per layer) did not have a massive effect on the final outcome. Specifically, increasing the number of layers or number of neurons per layer did not lead to large changes in model performance, though the training time increases rapidly as the number of connections increase. Code for, and results of, the final trained model are presented.
  * While initial models were very good at identifying control samples, the overall percentage of correct predictions was at most 70%. Closer investigation demonstrated a large number of the incorrect predictions stemmed from samples which were predicted as closely related MeSH terms or were samples which actually did appear in both groups. The model test function was expanded to allow for collection of the top X predictions for each sample, defaulting to the top 3 predictions. This led to a greater than 90% prediction accuracy.
    * The output file ’ DNN_output_epoch_500.csv’ demonstrates the potential output of the model and examples of a top 3 prediction. As expected, heavily related terms with overlapping sample groupings lead to a high level of incorrect first predictions but correct top 3 predictions.
    * Interestingly, [paranasal sinus diseases], [rhinitis], [allergic rhinitis], and [sinusitis] frequently demonstrated incorrect prediction even within the top 3 groups, predicting as some combination of [hypersensitivity, immediate], [lung disease, obstructive], and [asthma]. While it is not surprising  that there is poor predictivity between these groups due to the near complete overlap of sinus disease sample groups, it is surprising that none of the disease designations developed a strong enough signature to even appear in the top 3 predictions. It would be interesting to investigate this further due to the heavy overlap between sinus diseases and more general hypersensitivity/obstructive diseases/asthma. Maybe all of these groups should be collapsed into a pan-obstructive/hypersensitive signature group.

<img width="914" height="343" alt="image" src="https://github.com/user-attachments/assets/3e6507e4-0f97-422e-aaef-afa347c54987" />

8.	Random Forest:
  * A random forest (RF) model was grown using the R library ranger which offers a highly optimized implementation of the random forest algorithm.
    * https://cran.r-project.org/web/packages/ranger/index.html
  * The RF model was grown with the following settings:
    * 2000 trees
    * Maximum depth of 20 splits
    * Number of factors assessed at each split of sqrt(gene number - 1)
    * Case weights set to 1 – (number of samples in group/total sample number)
    * Permutation importance
  * Permutation importance calculates the importance of each individual factor (gene) by randomly permutation that factor and assessing the change in model performance from that change. This provides a measure of how important an individual factor may be in predicting the sample class.
    * Overall accuracy was quite low with similar performance to the DNN. When considering class probabilities rather than class labels, 99.1% accuracy was achieved when considering predictions within the top 3 highest predictive probabilities.
<img width="781" height="79" alt="image" src="https://github.com/user-attachments/assets/8af1b94e-5cf7-4f0b-8f77-7a64f408d54d" />
