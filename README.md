# Gottwein_CRISPR_Variant_Calling
Script for characterizing and quantiftying variants based on data returned from the Massachusetts General Hospital CCIB DNA Core.

# Requirements
You will need a version of Python supporting typehinting
pandas, biopython, crispresso2 (CLI tool)

# Usage
Generate a run table following the example provided

Perform run by passing run table path and the maximum number of codons that should be allowed out of frame before the gene is considered frameshifted (only relevant if you want to ignore when the frame is restored by a downstream indel very late in the read)
i.e. 
```
python ./crispr_characterizer.py run_table.csv 50
```

You will need to do this in an environment where CRISPResso2 is installed and accessible via the CRISPResso command.
