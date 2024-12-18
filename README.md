# Korean Blood Transcriptomic Clock

![image](https://github.com/user-attachments/assets/1621a068-c7cb-4d9e-a898-6b73d3b16093)


# Pre-processing
1. normalize_expression_data_train.R
2. filter_out_samples_by_age_group.py
3. generate_feature_table_gene_expression_per_age.py
4. stratify_split_train_and_test_feature_table.py
5. (fix_10s_80s_in_train_20s_in_test.py)
6. remove_median_zero_expression_train.py
7. submit_multiple_run_correlation_age.py

# Training mRNA Clock
1. standardize_testing_data.py
2. **make_transcriptomic_clock.py**

# Testing mRNA Clock
1. get_geoMeans_external_data.py
2. normalize_external_expression_data.R
3. generate_feature_table_gene_expression_per_age.py
4. standardize_testing_data.py
5. **KoreanBloodClock.py**

# Testing Other mRNA Clocks
- PetersClock.py
- RenClock.py

# Source 
- utility.py
- transcriptomic_clock.py
- run_correlation_age_gene_expression.py
- draw_enrichment_plot.py
- compare_clinical_values.py
- missing.py

# Requirements
conda environment
```
conda env create --file environment.yaml
```

R packages 
- DESeq2
- stringr
- argparse
