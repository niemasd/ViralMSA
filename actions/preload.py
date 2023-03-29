#! /usr/bin/env python3
import os
import sys
sys.path.append('./')

import ViralMSA

for key, value in ViralMSA.REFS.items():
    ViralMSA.download_ref_genome(value, 'website/assets/indexes/' + key, 'website/assets/indexes/' + key + '/' + key + '.fas', 'email@address.com')
    ViralMSA.build_index_minimap2('website/assets/indexes/' + key + '/' + key + '.fas', 1)