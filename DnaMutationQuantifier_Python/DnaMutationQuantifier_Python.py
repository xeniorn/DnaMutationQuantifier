from __future__ import annotations

# type hinting
import typing
from typing import List


# 0 - mostly silent, 1 - some stuff, 2 - debugging
global_verbosity = 2

def log(message, message_verbosity : int = -1):
  if message_verbosity <= global_verbosity:
    print(message)

from Bio.SeqIO.FastaIO import FastaTwoLineParser
#@title Shared object setup

# reading abi files
from Bio import SeqIO
#plotting
import matplotlib.pyplot as plt
# general
import math
from collections import defaultdict
#dna
from Bio import Seq

class Dna:
  @staticmethod
  def reverse_complement(sequence : str) -> str:
    return Seq.reverse_complement(sequence)

class SequencingLocus:
    def __init__(self, a, c, g, t):
        val_sum = a + c + g + t
        self.a = 0
        self.c = 0
        self.g = 0
        self.t = 0
        if val_sum == 0:
          return
        self.a = a / val_sum
        self.c = c / val_sum
        self.g = g / val_sum
        self.t = t / val_sum

    def __getitem__(self, key):
        if key == "A":
            return self.a
        elif key == "C":
            return self.c
        elif key == "G":
            return self.g
        elif key == "T":
            return self.t

    @staticmethod
    def get_complement(source_locus : SequencingLocus) -> SequencingLocus:
      return SequencingLocus(
          a = source_locus.t,
          c = source_locus.g,
          g = source_locus.c,
          t = source_locus.a
      )

class SangerTraceData:

  def __init__(self, a_data, c_data, g_data, t_data):
    self.data = {
        "A" : a_data,
        "C" : c_data,
        "G" : g_data,
        "T" : t_data,
    }
    self.a_data = a_data
    self.c_data = c_data
    self.g_data = g_data
    self.t_data = t_data

  def __getitem__(self, key):
    return self.data[key]


  def get_subset(self, start : int, end : int) -> SangerTraceData :
    s = slice(start, end-1)
    return SangerTraceData(self.data["A"][s], self.data["C"][s], self.data["G"][s], self.data["T"][s])

class NormalizedTrace:

    @staticmethod
    def from_abi_file(abi_file : str) -> NormalizedTrace:

        obj = NormalizedTrace()

        obj.is_reversed = False
        obj.file = abi_file

        record = SeqIO.read(abi_file, "abi")

        list(record.annotations.keys())
        list(record.annotations["abif_raw"].keys())

        order = record.annotations['abif_raw']['FWO_1']

        base_to_data = {
        order[0:1].decode(): "DATA9",
        order[1:2].decode(): "DATA10",
        order[2:3].decode(): "DATA11",
        order[3:4].decode(): "DATA12"
        }

        channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
        trace = defaultdict(list)
        for c in channels:
            trace[c] = record.annotations["abif_raw"][c]

        indices = record.annotations["abif_raw"]["PLOC2"]
        primary_basecall = record.annotations["abif_raw"]["PBAS2"]
        quality = list(record.annotations["abif_raw"]["PCON2"])

        max_qual = max(quality)

        obj.quality_absolute_max = max_qual

        obj.peak_indices = indices
        obj.sequence = primary_basecall.decode()
        obj.quality = [q/max_qual for q in quality]

        seqlen = len(obj.sequence)

        obj.raw_data = SangerTraceData(
            list(trace[base_to_data["A"]]),
            list(trace[base_to_data["C"]]),
            list(trace[base_to_data["G"]]),
            list(trace[base_to_data["T"]]))

        complex_trace = []

        for index in range(seqlen):
            peak_index = indices[index]
            vA = obj.raw_data["A"][peak_index]
            vC = obj.raw_data["C"][peak_index]
            vG = obj.raw_data["G"][peak_index]
            vT = obj.raw_data["T"][peak_index]
            seqloc = SequencingLocus(vA,vC,vG,vT)
            complex_trace.append(seqloc)

        obj.complex_trace = complex_trace

        log(f"Imported trace with {len(complex_trace)} points",2)

        return obj

    def __getitem__(self, key):
        return self.complex_trace[key]

    def calc_match(self, search_seq, template_start_index):
        seqlen = len(search_seq)

        if template_start_index + len(search_seq) > len(self.sequence):
            return 0

        score=0
        for probe_index in range(seqlen):
            template_index = template_start_index + probe_index
            score+=self[template_index][search_seq[probe_index]]
        return score

    def calc_matches(self, search_seq):

        seqlen = len(self.complex_trace)
        scores=[self.calc_match(search_seq, index) for index in range(seqlen)]

        return scores

    def get_top_match(self, search_seq : str, min_score_pc : float = 0.5) -> int:

        matches = self.calc_matches(search_seq)

        top_score = max(matches)
        log(f"top_score is {top_score}",3)
        if top_score < min_score_pc:
          return -1
        else:
          return matches.index(top_score)

    def get_value(self, nucleotide : str, nt_index : int) -> int:
      data_index = self.peak_indices[nt_index]
      return self.raw_data[nucleotide][data_index]

    @staticmethod
    def reverse(source_nt : NormalizedTrace) -> NormalizedTrace:
      result_nt = NormalizedTrace()

      result_nt.is_reversed = True
      result_nt.file = source_nt.file
      result_nt.quality_absolute_max = source_nt.quality_absolute_max
      result_nt.quality = list(reversed(source_nt.quality))

      result_nt.sequence = Dna.reverse_complement(source_nt.sequence)

      min_index = 0
      max_index = len(source_nt.raw_data.a_data) - 1
      result_nt.peak_indices = list(reversed([min_index + (max_index-ind) for ind in source_nt.peak_indices]))

      result_nt.raw_data = SangerTraceData(
          a_data = list(reversed(source_nt.raw_data.t_data)),
          c_data = list(reversed(source_nt.raw_data.g_data)),
          g_data = list(reversed(source_nt.raw_data.c_data)),
          t_data = list(reversed(source_nt.raw_data.a_data)))

      result_nt.complex_trace = list(reversed([SequencingLocus.get_complement(locus) for locus in source_nt.complex_trace]))

      return result_nt

class RatioToEditingConverter:
  def __init__(self, conversion_matrix):

    self.ratios = []
    self.efficiencies = []

    for entry in conversion_matrix:
      self.ratios.append(entry[0])
      self.efficiencies.append(entry[1])

  def convert(self, ratio : float) -> float:

    has_exact = self.ratios.count(ratio) > 0

    if has_exact:
      index = self.ratios.index(ratio)
      return self.efficiencies[index]

    if ratio < min(self.ratios):
      return -1

    if ratio > max(self.ratios):
      return -1

    # interpolate between the closest two values
    si = max([ index for index, r in enumerate(self.ratios) if r < ratio ])
    bi = min([ index for index, r in enumerate(self.ratios) if r > ratio ] )

    sr = self.ratios[si]
    br = self.ratios[bi]

    se = self.efficiencies[si]
    be = self.efficiencies[bi]

    final_eff = se + (be-se)*(br-ratio)/(br-sr)

    return final_eff



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

  def __init__(self, file : str = "", site : str = "", data : SangerTraceData = None, is_reverse : bool = False):
    self.file = file
    self.site = site
    self.data = data
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
  def extract_adenosine_substitution_values(nt : NormalizedTrace, site : str, converter : RatioToEditingConverter, is_reverse : bool ) -> List[TraceQuantRecord]:

    index = nt.get_top_match(site)
    adenosine_loc = [i for i, x in enumerate (site) if x == "A"]
    records = []

    for loc in adenosine_loc:
      abi_nt_index = index + loc
      value_g = nt[abi_nt_index]["G"]
      value_a = nt[abi_nt_index]["A"]
      log(f"G: {value_g} A: {value_a}",3)
      #value_g = nt.get_value("G", abi_nt_index)
      #value_a = nt.get_value("A", abi_nt_index)

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

    sanger_data = nt.raw_data.get_subset(start_subtrace_index, end_subtrace_index)

    res = SangerTraceDataRecord(nt.file, site, sanger_data, is_reverse)
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

    target = experiment.grna_with_target_setup.target_site
    off_sites = experiment.grna_with_target_setup.off_sites
    conv = experiment.ratio_to_editing_converter

    log(f"Got {len(experiment.controls)} controls and {len(experiment.targets)} targets.",1)
    log(f"Sites are: ", 2)
    log(target, 2)
    for off_site in off_sites:
      log(off_site, 2)

    results = []
    subtrace_records = []

    for ntc in experiment.controls:
      log(f"Working on control file {ntc.file}...",1)
      dir = self.is_it_fwd_or_reverse(ntc)
      log(f"Dir is {dir}", 1)
      if (dir == "fwd"):
        results.extend(self.get_results_for_existing_sites(ntc))
        subtrace_records.extend(self.get_subtraces_for_existing_sites(ntc))
      elif (dir == "rev"):
        ntc_rev = NormalizedTrace.reverse(ntc)
        results.extend(self.get_results_for_existing_sites(ntc_rev, True))
        subtrace_records.extend(self.get_subtraces_for_existing_sites(ntc_rev, True))
      else:
          log(f"Couldn't find unambiguous directionality for trace {ntc.file}", 1)

    for ntt in experiment.targets:
      log(f"Working on target file {ntt.file}...",1)
      dir = self.is_it_fwd_or_reverse(ntt)
      log(f"Dir is {dir}", 1)
      if (dir == "fwd"):
        results.extend(self.get_results_for_existing_sites(ntt))
        subtrace_records.extend(self.get_subtraces_for_existing_sites(ntt))
      elif (dir == "rev"):
        ntt_rev = NormalizedTrace.reverse(ntt)
        results.extend(self.get_results_for_existing_sites(ntt_rev, True))
        subtrace_records.extend(self.get_subtraces_for_existing_sites(ntt_rev, True))
      else:
          log(f"Couldn't find unambiguous directionality for trace {ntt.file}", 1)

    self.quant_results = results
    self.subtrace_results = subtrace_records

  @staticmethod
  def generate_standard_plots(data : SangerTraceData, ax : plt.Axes):
    ax.plot(data["A"], color="red", linewidth = 0.5)
    ax.plot(data["C"], color="blue", linewidth = 0.5)
    ax.plot(data["G"], color="black", linewidth = 0.5)
    ax.plot(data["T"], color="green", linewidth = 0.5)

  @staticmethod
  def plot_subtrace(rec : SangerTraceDataRecord, row : int, column : int, fig : plt.Figure, gs : plt.GridSpec):

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

  def generate_raw_plots(self):

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
        ExperimentAnalysis.plot_subtrace(rec, ifile, isite, fig, gs)

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

# mapping of G-to-A peak ratios to "editing efficiency", reformatted from the table in the .docx file that was provided
conversion_matrix = [
[0.0000,0.0000],
[0.0083,0.0080],
[0.0167,0.0160],
[0.0250,0.0240],
[0.0333,0.0320],
[0.0417,0.0400],
[0.0500,0.0480],
[0.0583,0.0550],
[0.0667,0.0630],
[0.0750,0.0700],
[0.0833,0.0770],
[0.0917,0.0840],
[0.1000,0.0910],
[0.1083,0.0980],
[0.1167,0.1040],
[0.1250,0.1110],
[0.1333,0.1180],
[0.1417,0.1240],
[0.1500,0.1300],
[0.1583,0.1370],
[0.1667,0.1430],
[0.1750,0.1490],
[0.1833,0.1550],
[0.1917,0.1610],
[0.2000,0.1670],
[0.2083,0.1720],
[0.2167,0.1780],
[0.2250,0.1840],
[0.2333,0.1890],
[0.2417,0.1950],
[0.2500,0.2000],
[0.2583,0.2050],
[0.2667,0.2110],
[0.2750,0.2160],
[0.2833,0.2210],
[0.2917,0.2260],
[0.3000,0.2310],
[0.3083,0.2360],
[0.3167,0.2410],
[0.3250,0.2450],
[0.3333,0.2500],
[0.3433,0.2550],
[0.3500,0.2590],
[0.3600,0.2640],
[0.3667,0.2680],
[0.3767,0.2730],
[0.3833,0.2770],
[0.3933,0.2810],
[0.4000,0.2860],
[0.4100,0.2900],
[0.4167,0.2940],
[0.4267,0.2980],
[0.4333,0.3020],
[0.4433,0.3060],
[0.4500,0.3100],
[0.4600,0.3140],
[0.4667,0.3180],
[0.4767,0.3220],
[0.4833,0.3260],
[0.4933,0.3300],
[0.5000,0.3330],
[0.5100,0.3370],
[0.5167,0.3410],
[0.5267,0.3440],
[0.5333,0.3480],
[0.5433,0.3510],
[0.5500,0.3550],
[0.5600,0.3580],
[0.5667,0.3620],
[0.5767,0.3650],
[0.5833,0.3680],
[0.5933,0.3720],
[0.6000,0.3750],
[0.6100,0.3780],
[0.6167,0.3810],
[0.6267,0.3850],
[0.6333,0.3880],
[0.6433,0.3910],
[0.6500,0.3940],
[0.6600,0.3970],
[0.6667,0.4000],
[0.6750,0.4030],
[0.6833,0.4060],
[0.6917,0.4090],
[0.7000,0.4120],
[0.7083,0.4150],
[0.7167,0.4170],
[0.7250,0.4200],
[0.7333,0.4230],
[0.7417,0.4260],
[0.7500,0.4290],
[0.7583,0.4310],
[0.7667,0.4340],
[0.7750,0.4370],
[0.7833,0.4390],
[0.7917,0.4420],
[0.8000,0.4440],
[0.8083,0.4470],
[0.8167,0.4500],
[0.8250,0.4520],
[0.8333,0.4550],
[0.8417,0.4570],
[0.8500,0.4590],
[0.8583,0.4620],
[0.8667,0.4640],
[0.8750,0.4670],
[0.8833,0.4690],
[0.8917,0.4710],
[0.9000,0.4740],
[0.9083,0.4760],
[0.9167,0.4780],
[0.9250,0.4810],
[0.9333,0.4830],
[0.9417,0.4850],
[0.9500,0.4870],
[0.9583,0.4890],
[0.9667,0.4920],
[0.9750,0.4940],
[0.9833,0.4960],
[0.9917,0.4980],
[1.0000,0.5000],
[1.0083,0.5020],
[1.0167,0.5040],
[1.0250,0.5060],
[1.0333,0.5080],
[1.0417,0.5100],
[1.0500,0.5120],
[1.0583,0.5140],
[1.0667,0.5160],
[1.0750,0.5180],
[1.0833,0.5200],
[1.0917,0.5220],
[1.1000,0.5240],
[1.1083,0.5260],
[1.1167,0.5280],
[1.1250,0.5290],
[1.1333,0.5310],
[1.1417,0.5330],
[1.1500,0.5350],
[1.1583,0.5370],
[1.1667,0.5380],
[1.1750,0.5400],
[1.1833,0.5420],
[1.1917,0.5440],
[1.2000,0.5450],
[1.2083,0.5470],
[1.2167,0.5490],
[1.2250,0.5510],
[1.2333,0.5520],
[1.2417,0.5540],
[1.2500,0.5560],
[1.2583,0.5570],
[1.2667,0.5590],
[1.2750,0.5600],
[1.2833,0.5620],
[1.2917,0.5640],
[1.3000,0.5650],
[1.3083,0.5670],
[1.3167,0.5680],
[1.3250,0.5700],
[1.3333,0.5710],
[1.3400,0.5730],
[1.3500,0.5740],
[1.3567,0.5760],
[1.3667,0.5770],
[1.3733,0.5790],
[1.3833,0.5800],
[1.3900,0.5820],
[1.4000,0.5830],
[1.4067,0.5850],
[1.4167,0.5860],
[1.4233,0.5880],
[1.4333,0.5890],
[1.4400,0.5900],
[1.4500,0.5920],
[1.4567,0.5930],
[1.4667,0.5950],
[1.4733,0.5960],
[1.4833,0.5970],
[1.4900,0.5990],
[1.5000,0.6000],
[1.5067,0.6010],
[1.5167,0.6030],
[1.5233,0.6040],
[1.5333,0.6050],
[1.5400,0.6070],
[1.5500,0.6080],
[1.5567,0.6090],
[1.5667,0.6100],
[1.5733,0.6120],
[1.5833,0.6130],
[1.5900,0.6140],
[1.6000,0.6150],
[1.6067,0.6170],
[1.6167,0.6180],
[1.6233,0.6190],
[1.6333,0.6200],
[1.6400,0.6210],
[1.6500,0.6230],
[1.6567,0.6240],
[1.6667,0.6250],
[1.6750,0.6260],
[1.6833,0.6270],
[1.6917,0.6280],
[1.7000,0.6300],
[1.7083,0.6310],
[1.7167,0.6320],
[1.7250,0.6330],
[1.7333,0.6340],
[1.7417,0.6350],
[1.7500,0.6360],
[1.7583,0.6370],
[1.7667,0.6390],
[1.7750,0.6400],
[1.7833,0.6410],
[1.7917,0.6420],
[1.8000,0.6430],
[1.8083,0.6440],
[1.8167,0.6450],
[1.8250,0.6460],
[1.8333,0.6470],
[1.8417,0.6480],
[1.8500,0.6490],
[1.8583,0.6500],
[1.8667,0.6510],
[1.8750,0.6520],
[1.8833,0.6530],
[1.8917,0.6540],
[1.9000,0.6550],
[1.9083,0.6560],
[1.9167,0.6570],
[1.9250,0.6580],
[1.9333,0.6590],
[1.9417,0.6600],
[1.9500,0.6610],
[1.9583,0.6620],
[1.9667,0.6630],
[1.9750,0.6640],
[1.9833,0.6650],
[1.9917,0.6660],
[2.0000,0.6670]
]

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

def fix_path(path: str):
  return path.replace("\\","/").replace('"', '')

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

controls = [NormalizedTrace.from_abi_file(fix_path(x)) for x in control_files]
targets = [NormalizedTrace.from_abi_file(fix_path(x)) for x in target_files]

grna_setup = GrnaWithTargetSetup(grna_sequence, target_sequence, list(off_sites.values()))
conv = RatioToEditingConverter(conversion_matrix)
experiment = QuantificationExperimentSetup(grna_setup, controls, targets, conv)

global_verbosity=1

anal = ExperimentAnalysis(experiment)

log("Processing...")
anal.process()
log("Generating plots...")
anal.generate_raw_plots()
log("Exporting quant...")
anal.export_quantification()