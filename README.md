# Processing & Filtering
Much of this project was done on the cloud via Terra.bio (Firecloud v2) 
workflows or Google VMs. As such, if you need access to the Terra workspaces, 
please email twood@broadinstitute.org.


User Warning: Many of the bash scripts used to train and
evaluate the models use nohup commands - if your CPU is not able to tolerate
the amount of jobs, consider modifying your local version to serially run the models.

# Generation of the Gene Sample Matrix (GSM, updated 24 Sep 2025)
Four python scripts in subdirectory src_python/gsm_generation_scripts translate somatic 
variants in common formats into a GSM with samples listed horizontally and DLBCL driver 
variants  vertically. The GSM can then be input to the downloadable DLBclass app 
(https://data.broadinstitute.org/dlbclass) which assigns each sample to one of 
five possible classes. 

1) maf2gsm.py: converts somatic mutations from Mutation Annotation Format (MAF) format 
(https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) to a partial GSM specific to 
mutations.  
2) sv2gsm.py: converts translocations from a tsv format similar to bedpe 
(https://bedtools.readthedocs.io/en/latest/content/general-usage.html). Optionally sv2gsm.py 
can also load translocations from a simple GSM-like format to create a combined 
translocation-specific GSM. 
3) seg2gsm.py: converts Copy Number Variants from seg format 
(https://igv.org/doc/desktop/#FileFormats/DataTracks/) to a CNV-specific GSM. 
4) combine2gsm.py: collects output from maf2gsm, sv2gsm, and seg2gsm to create a unified 
GSM including rows for all 163 DLBCL driver variants used by DLBclass. 

# Terra Workflow Links

#### Main pipeline (Mutect, Mutsig2CV, Strelka, & others)

#### CGA HG38A WES Pipeline, split into pair sets :

NCI_A: https://app.terra.bio/#workspaces/shipp-dfci/hg38_tumorOnly_test_workspace/job_history/b0649a72-3c73-412f-b895-3d7779d43b88

NCI_B: https://app.terra.bio/#workspaces/shipp-dfci/hg38_tumorOnly_test_workspace/job_history/35b939df-447e-469f-8c62-e6829b366ab4

NCI_C: https://app.terra.bio/#workspaces/shipp-dfci/hg38_tumorOnly_test_workspace/job_history/a9d0751e-4be1-4d42-9b32-ada9637874a4

NCI_D: https://app.terra.bio/#workspaces/shipp-dfci/hg38_tumorOnly_test_workspace/job_history/8e504071-f828-4360-be1f-2aed5d3fca6d

NCI_E: https://app.terra.bio/#workspaces/shipp-dfci/hg38_tumorOnly_test_workspace/job_history/ddf231d9-5c95-48a6-afcc-673ad5b48649

#### Standalone BLAT filter and beyond (call-cached jobs, in order)

https://app.terra.bio/#workspaces/shipp-dfci/DLBCL_Staudt_TumorOnly_2021_v2/job_history

# Classifier Preprocessing Reproducibility Steps (optional, files already included)

1) Remap the labels from the consensus nmf job above: src_python/remap_labels.py
2) Compute q-values per gene: src_python/fisher_5x2_parallel.py
3) Generate baseline probabilities: src_python/generate_baseline_probabilities.py
4) Create gene footprint table: src_python/calculate_driver_footprint.py

# Before training Classifier:

#### Create an environment

conda create --name Classifier

conda activate Classifier

##### Install Pre-requisites

conda install pytorch torchvision -c pytorch

conda install pandas

conda install matplotlib

conda install scikit-learn

# Model Reproducibility steps

1) Run model training bash scripts (warning: this will launch many jobs, do not launch at all once)

    * run_all_experiments_step1.sh
    
    * run_all_experiments_step2A.sh
    
    * run_all_experiments_step2B.sh
    
    * run_all_experiments_step2C.sh
    
    * run_all_experiments_step2T.sh
    
    * run_sens_spec_experiments.sh

2) Evaluate all trained models: src_python/evaluate_validation_ensembles.py
3) Combine training history: src_python/combine_model_training_history.py

# Plotting Results

Most plots are generated via R

1) Step 1: src_R/plot_step1.R
2) Step 2A: src_R/plot_step2A.R
3) Step 2B: src_R/plot_step2B.R
4) Step 2C: src_R/plot_step2C.R
5) Step 2T: src_R/plot_step2T.R
6) Sens/Spec experiments: src_R/plot_sensitivity_specificity_experiment.R
7) Model training history: src_R/plot_training_history.R

