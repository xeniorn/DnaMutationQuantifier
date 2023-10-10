from dna_mutation_quantifier.core import log
from dna_mutation_quantifier.run import default_process_experiment, init_experiment


target_sequence = 'TGAAGCGCACGTAGGAGCAG' #@param {type:"string"}
#@markdown
recruitment_cluster_1 = '' #@param {type:"string"}
recruitment_cluster_2 = '' #@param {type:"string"}
recruitment_cluster_3 = '' #@param {type:"string"}
recruitment_cluster_4 = '' #@param {type:"string"}

#@markdown (not used for now):
grna_sequence = 'XXXXX' #@param {type:"string"}

off_sites = {
    1: recruitment_cluster_1,
    2: recruitment_cluster_2,
    3: recruitment_cluster_3,
    4: recruitment_cluster_4,
  }

min_cluster_length=4

off_sites = {key: seq.upper() for key, seq in off_sites.items() if len(seq) >= min_cluster_length}

print(off_sites)

control_files = [
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_mb18525.ab1"',
]

target_files = [
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_V2_IFN_mb18525.ab1"',
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_V2_mb18525.ab1"',
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_V3_IFN_mb18525.ab1"',
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_V3_mb18525.ab1"',
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_V1_IFN_mb18525.ab1"',
r'"C:\temp\dnaquant\ADNP-PTC719(TAG)_B4_V1_mb18525.ab1"',
]

#target_files =[]

mode="precise"
mode="rapid"

use_calibration_curve = False

show_peaks = False

anal = init_experiment(control_files, target_files, grna_sequence, target_sequence, off_sites, mode, use_calibration_curve)

log("Processing...")
anal.process()
log("Generating plots...")
anal.generate_raw_plots(show_peaks)
log("Exporting quant...")
anal.export_quantification()

#default_process_experiment(anal)