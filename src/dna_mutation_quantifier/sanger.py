from __future__ import annotations

from dna_mutation_quantifier.core import log, Constants

# type hinting
import typing
from typing import List

# reading abi files
from Bio import SeqIO

# fasta parsing
from Bio import Seq
from Bio.SeqIO.FastaIO import FastaTwoLineParser

# general
import math
import numpy as np
import scipy as sp
from collections import defaultdict

class MathHelpers:
    def auto_fit_gaussian(data_slice: List[float]):
                    
                    def gaussian(x: float, A: float, mu: float, sigma: float, offset: float) -> float:
                        return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + offset
                    
                    x = np.linspace(0, len(data_slice)-1, len(data_slice))
                    y = np.array(data_slice)

                    from scipy.optimize import curve_fit

                    init_guess = np.array([max(data_slice), len(data_slice)/2, len(data_slice)/3, 0])
                    try:
                        params, covariance = curve_fit(gaussian, x, y, init_guess)
                        A_fit, mu_fit, sigma_fit, offset_fit = params
                        return (A_fit, mu_fit, sigma_fit, offset_fit, covariance)
                    except:
                        return None


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


  def get_subset(self, start : int, end : int) -> SangerTraceData:
    s = slice(start, end-1)
    return SangerTraceData(self.data["A"][s], self.data["C"][s], self.data["G"][s], self.data["T"][s])

class NormalizedTrace:    

    def __init__(self, 
                 raw_data: SangerTraceData = None,
                 raw_peaks: typing.Dict[str, List[int]] = None):
        self.raw_data = raw_data
        self.raw_peaks = raw_peaks

    def get_subset(self, start_subtrace_index: int, end_subtrace_index: int) -> typing.Tuple[SangerTraceData, typing.Dict[str, List[int]]]:
        data = self.raw_data.get_subset(start_subtrace_index, end_subtrace_index)
        def get_peaks_in_subtrace(nucl_peaks: List[int]) -> List[int]:
            return [peak - start_subtrace_index for peak in nucl_peaks if start_subtrace_index <= peak <= end_subtrace_index]
        peaks = {nucl: get_peaks_in_subtrace(nucl_peaks) for nucl, nucl_peaks in self.raw_peaks.items()}
        return (data, peaks)
    
    def refine_complex_trace(start_nt_index: int, end_nt_index: int) -> List[SangerTraceData]:
        return []

    @staticmethod
    def from_abi_file(abi_file : str, mode: str = "rapid") -> NormalizedTrace:
        """
        abi_file:
        path to .ab1 file from which the trace will be constructed
        
        mode:        
        "rapid" - get values for each nucleotide from the annotated peak
        "precise" - do peak detection for each nucleotide channel signal and quantify at peak
        """

        obj = NormalizedTrace()

        obj.is_reversed = False
        obj.file = abi_file

        record = SeqIO.read(abi_file, "abi")

        list(record.annotations.keys())
        list(record.annotations["abif_raw"].keys())

        order = record.annotations['abif_raw']['FWO_1']

        # maps ACGT to data fields
        base_to_data = {
        order[0:1].decode(): "DATA9",
        order[1:2].decode(): "DATA10",
        order[2:3].decode(): "DATA11",
        order[3:4].decode(): "DATA12"
        }
       

        channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
        
        trace = {c: record.annotations["abif_raw"][c] for c in channels}

        #trace = defaultdict(list)
        #for c in channels:
        #    trace[c] = record.annotations["abif_raw"][c]

        obj.peak_indices = record.annotations["abif_raw"]["PLOC2"]
        primary_basecall = record.annotations["abif_raw"]["PBAS2"]
        quality = list(record.annotations["abif_raw"]["PCON2"])

        obj.quality_absolute_max = max(quality)
        obj.quality = [q/obj.quality_absolute_max for q in quality]

        obj.sequence = primary_basecall.decode()
        
        obj.raw_data = SangerTraceData(
            list(trace[base_to_data["A"]]),
            list(trace[base_to_data["C"]]),
            list(trace[base_to_data["G"]]),
            list(trace[base_to_data["T"]]))
        
        obj.raw_data_max_index = len(obj.raw_data.a_data) - 1
                
        if mode == "rapid":
            obj.set_complex_trace_rapid()
        elif mode == "precise":
            obj.set_complex_trace_precise()
        else:
            raise Exception(f"Unsupported trace quantification mode: {mode}")
        
        log(f"Imported trace with {len(obj.complex_trace)} points",2)

        return obj
    
    def set_complex_trace_rapid(self):

        complex_trace = []
        seqlen = len(self.sequence)

        self.raw_peaks = {nucl: [] for nucl in Constants.nucleotides}

        for index in range(seqlen):
            peak_index = self.peak_indices[index]
            vA = self.raw_data["A"][peak_index]
            vC = self.raw_data["C"][peak_index]
            vG = self.raw_data["G"][peak_index]
            vT = self.raw_data["T"][peak_index]
            seqloc = SequencingLocus(vA,vC,vG,vT)
            complex_trace.append(seqloc)

            for nucl in Constants.nucleotides:
                self.raw_peaks[nucl].append(peak_index)            

        self.complex_trace = complex_trace

    def set_complex_trace_precise(self):        
        
        narrowing_factor = Constants.peak_search_narrowing_factor

        complex_trace = []
        seqlen = len(self.sequence)

        self.raw_peaks = {nucl: [] for nucl in Constants.nucleotides}

        for index in range(seqlen):           
            
            if index == 150:
                print("150!")

            peaks = self.peak_indices

            peak_curr = peaks[index]
            
            left_bound_raw = (peaks[index-1] + peak_curr)/2 if index > 0 else 0
            right_bound_raw = (peaks[index+1] + peak_curr)/2 if index < seqlen-1 else self.raw_data_max_index

            left_bound = math.trunc(peak_curr - (peak_curr - left_bound_raw) * narrowing_factor)
            right_bound = math.trunc(peak_curr + (right_bound_raw - peak_curr) * narrowing_factor)

            # def find_max_index_gauss_BROKEN(nucleotide: str) -> int:
            #     data_slice = self.raw_data[nucleotide][left_bound:right_bound+1]

            #     fit_res = MathHelpers.auto_fit_gaussian(data_slice)
            #     if fit_res is not None:
            #         (A_fit, mu_fit, _sigma_fit, _offset_fit, _covariance) = fit_res
            #         max_index = round(left_bound + mu_fit)                
            #         return max_index if left_bound<=max_index<=right_bound else peak_curr
            #     else:
            #         return peak_curr
            
            def find_max_index_simple(nucleotide: str) -> int:
                data_slice = self.raw_data[nucleotide][left_bound:right_bound+1]
                max_val = max(data_slice)
                max_index = data_slice.index(max_val)                
                return max_index


            nuclPeaks = {nucl: find_max_index_simple(nucl) for nucl in Constants.nucleotides}
            nuclValues = {nucl: self.raw_data[nucl][nuclPeaks[nucl]] for nucl in Constants.nucleotides}
            
            seqloc = SequencingLocus(nuclValues["A"],nuclValues["C"],nuclValues["G"],nuclValues["T"])
            complex_trace.append(seqloc)

            for nucl in Constants.nucleotides:
                self.raw_peaks[nucl].append(nuclPeaks[nucl]) 
        
        self.complex_trace = complex_trace

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
      max_index = source_nt.raw_data_max_index
      result_nt.peak_indices = list(reversed([min_index + (max_index-ind) for ind in source_nt.peak_indices]))

      result_nt.raw_data = SangerTraceData(
          a_data = list(reversed(source_nt.raw_data.t_data)),
          c_data = list(reversed(source_nt.raw_data.g_data)),
          g_data = list(reversed(source_nt.raw_data.c_data)),
          t_data = list(reversed(source_nt.raw_data.a_data)))
      
      result_nt.raw_peaks = {
          "A": list(reversed([min_index + (max_index-ind) for ind in source_nt.raw_peaks["T"]])),
          "C": list(reversed([min_index + (max_index-ind) for ind in source_nt.raw_peaks["G"]])),
          "G": list(reversed([min_index + (max_index-ind) for ind in source_nt.raw_peaks["C"]])),
          "T": list(reversed([min_index + (max_index-ind) for ind in source_nt.raw_peaks["A"]])),          
      }

      result_nt.complex_trace = list(reversed([SequencingLocus.get_complement(locus) for locus in source_nt.complex_trace]))
    
      result_nt.raw_peaks 
    
      return result_nt

