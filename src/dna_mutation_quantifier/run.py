from dna_mutation_quantifier.core import *
from dna_mutation_quantifier.calibration import RatioToEditingConverter, conversion_matrix
from dna_mutation_quantifier.sanger import NormalizedTrace
from dna_mutation_quantifier.analysis import ExperimentAnalysis, GrnaWithTargetSetup, QuantificationExperimentSetup

import typing

def fix_path(path: str):
  return path.replace("\\","/").replace('"', '')

def init_experiment(control_files: typing.List[str], 
                    target_files: typing.List[str], 
                    grna_sequence: str, 
                    target_sequence: str, 
                    off_sites: typing.Dict[int,str], 
                    mode: str = "precise",
                    use_calibration_curve: bool = False):
    
    controls = [NormalizedTrace.from_abi_file(fix_path(x), mode) for x in control_files]
    targets = [NormalizedTrace.from_abi_file(fix_path(x), mode) for x in target_files]

    grna_setup = GrnaWithTargetSetup(grna_sequence, target_sequence, list(off_sites.values()))
    if use_calibration_curve:
        conv = RatioToEditingConverter(conversion_matrix)
    else:
        conv = RatioToEditingConverter()
    experiment = QuantificationExperimentSetup(grna_setup, controls, targets, conv)

    anal = ExperimentAnalysis(experiment)

    return anal


def default_process_experiment(anal: ExperimentAnalysis):  
    log("Processing...")
    anal.process()
    log("Generating plots...")
    anal.generate_raw_plots()
    log("Exporting quant...")
    anal.export_quantification()