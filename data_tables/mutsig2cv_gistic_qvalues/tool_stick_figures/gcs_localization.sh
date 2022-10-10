source '/cromwell_root/gcs_transfer.sh'

timestamped_message 'Localization script execution started...'


# No reference disks mounted since not requested in workflow options.





# Localize files from source bucket 'broad-institute-gdac' to container parent directory '/cromwell_root/broad-institute-gdac/reference/UniProt'.
files_to_localize_3e444e5c002c20baefa601337a6c3dcc=(
  "shipp-dfci"   # project to use if requester pays
  "3" # max transfer attempts
  "/cromwell_root/broad-institute-gdac/reference/UniProt/" # container parent directory
  "gs://broad-institute-gdac/reference/UniProt/Uniprot_SWISS_human_Release_2017_10_25.shelf"
  "gs://broad-institute-gdac/reference/UniProt/UniProtKB_AC_to_ID_human_Release_2017_10_25.shelf"
)

localize_files "${files_to_localize_3e444e5c002c20baefa601337a6c3dcc[@]}"
       



# Localize files from source bucket 'fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5' to container parent directory '/cromwell_root'.
files_to_localize_dbf3d6398678ffce589e606e65d25cbe=(
  "shipp-dfci"   # project to use if requester pays
  "3" # max transfer attempts
  "/cromwell_root/" # container parent directory
  "gs://fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5/submissions/c1fa8780-6e0e-4e55-afef-a2e39cf9194e/mutation_mutsig2cv_hg19/38c35235-dd76-4d8a-bc4b-55309cd27785/call-tool_stick_figures/script"
)

localize_files "${files_to_localize_dbf3d6398678ffce589e606e65d25cbe[@]}"
       


# Localize files from source bucket 'fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5' to container parent directory '/cromwell_root/fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5/submissions/c1fa8780-6e0e-4e55-afef-a2e39cf9194e/mutation_mutsig2cv_hg19/38c35235-dd76-4d8a-bc4b-55309cd27785/call-tool_mutsig2cv_hg19'.
files_to_localize_862c64b3e17168418ae8b9b7be1eb074=(
  "shipp-dfci"   # project to use if requester pays
  "3" # max transfer attempts
  "/cromwell_root/fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5/submissions/c1fa8780-6e0e-4e55-afef-a2e39cf9194e/mutation_mutsig2cv_hg19/38c35235-dd76-4d8a-bc4b-55309cd27785/call-tool_mutsig2cv_hg19/" # container parent directory
  "gs://fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5/submissions/c1fa8780-6e0e-4e55-afef-a2e39cf9194e/mutation_mutsig2cv_hg19/38c35235-dd76-4d8a-bc4b-55309cd27785/call-tool_mutsig2cv_hg19/sig_genes.txt"
  "gs://fc-secure-359a6793-cdf6-4459-a54d-d3a7481811e5/submissions/c1fa8780-6e0e-4e55-afef-a2e39cf9194e/mutation_mutsig2cv_hg19/38c35235-dd76-4d8a-bc4b-55309cd27785/call-tool_mutsig2cv_hg19/Combined_699_nopde4dip.final_analysis_set.maf"
)

localize_files "${files_to_localize_862c64b3e17168418ae8b9b7be1eb074[@]}"
       

timestamped_message 'Localization script execution complete.'