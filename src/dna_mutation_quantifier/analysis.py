from typing import List

from dna_mutation_quantifier.core import *
from dna_mutation_quantifier.sanger import NormalizedTrace, SangerTraceData
from dna_mutation_quantifier.calibration import *

import math
import matplotlib.pyplot as plt

class GrnaWithTargetSetup:
  def __init__(self, grna_seq : str, target_site : str, off_sites : List[str]):
    self.grna_seq = grna_seq
    self.target_site = target_site
    self.off_sites = off_sites

class QuantificationExperimentSetup:
  def __init__(self, grna_with_target_setup : GrnaWithTargetSetup, controls : List[NormalizedTrace], targets : List[NormalizedTrace], ratio_to_editing_converter):
    self.grna_with_target_setup = grna_with_target_setup
    self.controls = controls
    self.targets = targets
    self.ratio_to_editing_converter = ratio_to_editing_converter

class TraceQuantRecord:

  def __init__(self, file : str = "", site : str = "", sub_site_index : int = -1, ratio : float = 0, substitution_factor : float = 0, is_reverse : bool = False):
    self.file = file
    self.site = site
    self.index = sub_site_index
    self.ratio = ratio
    self.substitution_factor = substitution_factor
    self.is_reverse = is_reverse

  @classmethod
  def header_as_tsv(cls) -> str:
    return "\t".join(["source_file",
                      "target_site",
                      "index_in_site",
                      "ratio",
                      "substitution_factor",
                      "is_reverse"])

  def as_tsv(self) -> str:
    return "\t".join([self.file,
                      self.site,
                      str(self.index),
                      str(round(self.ratio,3)),
                      str(round(self.substitution_factor,3)),
                      str(self.is_reverse)])

class SangerTraceDataRecord:

  def __init__(self, file : str = "", site : str = "", data : SangerTraceData = None, peaks : typing.Dict[str, int] = None, is_reverse : bool = False):
    self.file = file
    self.site = site
    self.data = data
    self.peaks = peaks
    self.is_reverse = is_reverse






class TraceAnalysis:

  @staticmethod
  def check_site(nt : NormalizedTrace, site : str) -> int:
    index = nt.get_top_match(site)
    if index >= 0:
      #print(f"Found site at {index}")
      return index
    else:
      print(f"Site not found")
      return -1

  @staticmethod
  def extract_adenosine_substitution_values(nt : NormalizedTrace, site : str, converter : RatioToEditingConverter, is_reverse : bool, mode: str = "rapid" ) -> List[TraceQuantRecord]:

    index = nt.get_top_match(site)
    adenosine_loc = [i for i, x in enumerate (site) if x == "A"]
    records = []

    if mode == "rapid":
      values = nt.complex_trace
    elif mode == "precise":
      start = index
      end = start + len(site) - 1
      values = nt.refine_complex_trace(start, end)


    for loc in adenosine_loc:
      abi_nt_index = index + loc
      value_g = values[abi_nt_index]["G"]
      value_a = values[abi_nt_index]["A"]
      log(f"G: {value_g} A: {value_a}",3)      

      if value_a > 0:
        ratio = value_g / value_a
        sub_val = converter.convert(ratio)
      else:
        if value_g > 0:
          ratio = math.inf
          sub_val = 1
        else:
          ratio = 0
          sub_val = 0

      adenosine_locus_base_one = loc + 1
      res = TraceQuantRecord(nt.file, site, adenosine_locus_base_one, ratio, sub_val, is_reverse)
      records.append(res)

    return records

  @staticmethod
  def extract_subtrace(nt : NormalizedTrace, site : str, is_reverse : bool = False) -> SangerTraceDataRecord :

    start_base_index = nt.get_top_match(site)
    end_base_index = start_base_index + len(site) - 1

    start_peak_trace_index = nt.peak_indices[start_base_index]
    end_peak_trace_index = nt.peak_indices[end_base_index]

    # this will break if it's only 1 base, should make it get the "local average width" around the selection, not necessarily _in_ the selection
    average_peak_width = (end_peak_trace_index - start_peak_trace_index) / (end_base_index - start_base_index)
    half_width = average_peak_width/2

    extend = half_width * 1.1

    start_subtrace_index = max(math.trunc(start_peak_trace_index - extend), 0)
    end_subtrace_index = min(math.trunc(end_peak_trace_index + extend + 1), len(nt.raw_data["A"]))

    log(f"ex: {extend}, spti: {start_peak_trace_index}, epti: {end_peak_trace_index}, ssi: {start_subtrace_index}, esi: {end_subtrace_index}", 3)

    sanger_data, peaks = nt.get_subset(start_subtrace_index, end_subtrace_index)

    res = SangerTraceDataRecord(nt.file, site, sanger_data, peaks, is_reverse)
    return res


  @staticmethod
  def get_plot_data(self, index : int, length : int):
    subtrace = ""
    return subtrace


class ExperimentAnalysis:

  def __init__(self, experiment : QuantificationExperimentSetup):
    self.experiment = experiment
    self.quant_results = []
    self.subtrace_results = []
    self.result_files = []
    self.result_file_base="res_" + self.experiment.grna_with_target_setup.grna_seq + "_" + self.experiment.grna_with_target_setup.target_site


  def is_it_fwd_or_reverse(self, nt : NormalizedTrace, min_rel_score : int = 0.75) -> str :

    fwd_score_sum = 0
    rev_score_sum = 0
    total_len = 0

    nt_fwd = nt
    nt_rev = NormalizedTrace.reverse(nt_fwd)

    target = self.experiment.grna_with_target_setup.target_site
    off_sites = self.experiment.grna_with_target_setup.off_sites

    for site in [target] + off_sites:
      fwd_score_sum+=max(nt_fwd.calc_matches(search_seq=site))
      rev_score_sum+=max(nt_rev.calc_matches(search_seq=site))
      total_len+=len(site)

    log(f"Sequence match scores are fwd:{str(round(fwd_score_sum / total_len, 2))}, rev:{str(round(rev_score_sum / total_len, 2))}",3)

    max_score = total_len
    max_rel_score = max(fwd_score_sum, rev_score_sum) / max_score
    if max_rel_score < min_rel_score or fwd_score_sum == rev_score_sum:
      return "ambiguous"

    if fwd_score_sum > rev_score_sum:
      return "fwd"
    else:
      return "rev"


  def get_subtraces_for_existing_sites(self, nt : NormalizedTrace, is_reverse = False) -> List[SangerTraceDataRecord] :

    target = self.experiment.grna_with_target_setup.target_site
    off_sites = self.experiment.grna_with_target_setup.off_sites
    conv = self.experiment.ratio_to_editing_converter

    results = []

    log(f"In {nt.file}, getting plot data for target site {target}...",2)
    if TraceAnalysis.check_site(nt, target) > 0:
      tempres = TraceAnalysis.extract_subtrace(nt, target, is_reverse)
      results.append(tempres)

    for off_site in off_sites:
      log(f"In {nt.file}, getting plot data for off-site {off_site}...",2)
      if TraceAnalysis.check_site(nt, off_site) > 0:
        tempres = TraceAnalysis.extract_subtrace(nt, off_site, is_reverse)
        results.append(tempres)

    return results


  def get_results_for_existing_sites(self, nt : NormalizedTrace, is_reverse = False) -> List[TraceQuantRecord] :

    target = self.experiment.grna_with_target_setup.target_site
    off_sites = self.experiment.grna_with_target_setup.off_sites
    conv = self.experiment.ratio_to_editing_converter

    results = []

    log(f"In {nt.file}, quantifying target site {target}...",2)
    if TraceAnalysis.check_site(nt, target) > 0:
      tempres = TraceAnalysis.extract_adenosine_substitution_values(nt, target, conv, is_reverse)
      results.extend(tempres)

    for off_site in off_sites:
      log(f"In {nt.file}, quantifying off-site {off_site}...",2)
      if TraceAnalysis.check_site(nt, off_site) > 0:
        tempres = TraceAnalysis.extract_adenosine_substitution_values(nt, off_site, conv, is_reverse)
        results.extend(tempres)

    return results

  def process(self):

    target = self.experiment.grna_with_target_setup.target_site
    off_sites = self.experiment.grna_with_target_setup.off_sites
    conv = self.experiment.ratio_to_editing_converter

    log(f"Got {len(self.experiment.controls)} controls and {len(self.experiment.targets)} targets.",1)
    log(f"Sites are: ", 2)
    log(target, 2)
    for off_site in off_sites:
      log(off_site, 2)

    results = []
    subtrace_records = []

    def process_set(trace_set: List[NormalizedTrace], set_name: str):
      for ntx in trace_set:
        log(f"Working on {set_name} file {ntx.file}...",1)        
        dir = self.is_it_fwd_or_reverse(ntx)
        log(f"Dir is {dir}", 1)
                
        if (dir == "fwd"):
          reverse = False                
        elif (dir == "rev"):
          reverse = True
          ntx = NormalizedTrace.reverse(ntx)                   
        else:
            log(f"Couldn't find unambiguous directionality for trace {ntx.file}", 1)
            return
        
        results.extend(self.get_results_for_existing_sites(ntx, is_reverse=reverse))
        subtrace_records.extend(self.get_subtraces_for_existing_sites(ntx, is_reverse=reverse))

    process_set(self.experiment.controls, "control")
    process_set(self.experiment.targets, "target")    

    self.quant_results = results
    self.subtrace_results = subtrace_records

  @staticmethod
  def generate_standard_plots(data : SangerTraceData, ax : plt.Axes):
    ax.plot(data["A"], color="red", linewidth = 0.5)
    ax.plot(data["C"], color="blue", linewidth = 0.5)
    ax.plot(data["G"], color="black", linewidth = 0.5)
    ax.plot(data["T"], color="green", linewidth = 0.5)

  @staticmethod
  def plot_subtrace(rec : SangerTraceDataRecord, row : int, column : int, fig : plt.Figure, gs : plt.GridSpec, show_peaks : bool):

    ax = fig.add_subplot(gs[row,column])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    #ax.set_title(f"{rec.site} :: {rec.file}", fontsize=5)
    if row == 0:
      ax.set_title(f"{rec.site}", fontsize=8)
    if column == 0:
      ax.yaxis.set_visible(True)
      ax.yaxis.set_ticks([])
      ax.set_ylabel(rec.file, fontsize=5, )
    ExperimentAnalysis.generate_standard_plots(rec.data, ax)

    if show_peaks:
      for nucl in Constants.nucleotides:
        if rec.peaks is not None and nucl in rec.peaks:
          for peak_index in rec.peaks[nucl]:
            line_offset = 20
            center = rec.data[nucl][peak_index]
            ax.vlines(peak_index, center - line_offset, center + line_offset, color="black")
            #ax.annotate(".", (peak_index, rec.data[nucl][peak_index]))


  def generate_raw_plots(self, show_peaks: bool = True):

    num_offsites = len(self.experiment.grna_with_target_setup.off_sites)
    num_samples = len(self.experiment.targets)
    num_controls = len(self.experiment.controls)

    all_sites = [self.experiment.grna_with_target_setup.target_site] + list(self.experiment.grna_with_target_setup.off_sites)
    all_files = [x.file for x in self.experiment.controls] + [x.file for x in self.experiment.targets]

    fig = plt.figure(0, figsize = (3 * len(all_sites), 2 * len(all_files)))

    gs = plt.GridSpec(num_samples + num_controls, 1 + len(all_sites))


    for ifile in range(len(all_files)):
      file = all_files[ifile]
      for isite in range(len(all_sites)):
        site = all_sites[isite]

        # there can be... ooonly one...
        rec = [x for x in self.subtrace_results if (x.file == file) & (x.site == site)][0]
        ExperimentAnalysis.plot_subtrace(rec, ifile, isite, fig, gs, show_peaks)

    plot_result_name1 = self.result_file_base + "_plot.svg"
    plot_result_name2 = self.result_file_base + "_plot.png"

    if self.result_files.count(plot_result_name1) == 0:
      self.result_files.append(plot_result_name1)
    if self.result_files.count(plot_result_name2) == 0:
      self.result_files.append(plot_result_name2)

    plt.savefig(plot_result_name1)
    plt.savefig(plot_result_name2, dpi=300)

    plt.show()

  def export_quantification(self):

    quant_result_file_name = self.result_file_base + "_quant.tsv"

    if self.result_files.count(quant_result_file_name) == 0:
      self.result_files.append(quant_result_file_name)

    f = open(quant_result_file_name, "w")

    control_filenames = [x.file for x in self.experiment.controls]


    header = "\t".join(["is_control", TraceQuantRecord.header_as_tsv()])
    f.write( header + '\n')
    for rec in self.quant_results:
      is_control = control_filenames.count(rec.file) > 0
      data_line="\t".join([str(is_control), rec.as_tsv()])
      f.write(data_line + '\n')

    f.close()

  def export_zip(self, output_filename : str = "") -> str:

    import zipfile

    if output_filename=="":
      zip_filename=self.result_file_base + ".zip"
    else:
      zip_filename=output_filename

    zf = zipfile.ZipFile(zip_filename, "w")

    for file in self.result_files:
      zf.write(file)

    zf.close()

    return zip_filename
  
  #@title Conversion Matrix



