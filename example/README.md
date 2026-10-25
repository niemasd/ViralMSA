This directory contains example files to help you verify that ViralMSA is set up properly.

* [`example_hiv.fas`](example_hiv.fas) is a FASTA file containing 100 HIV-1 whole genome sequences. You can run ViralMSA as follows:
    * `ViralMSA.py -e email@address.com -s example_hiv.fas -o output -r HIV-1`
* [`reference_hiv.fas`](reference_hiv.fas) is a FASTA file containing the HIV-1 reference genome sequence. You can run ViralMSA as follows:
    * `ViralMSA.py -s example_hiv.fas -o output -r reference_hiv.fas`
