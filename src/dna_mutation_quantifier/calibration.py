import typing

class RatioToEditingConverter:
  def __init__(self, conversion_matrix = None):

    if conversion_matrix is None:
      self.calibrated = False
      return
    else:
      self.calibrated = True

    self.ratios = []
    self.efficiencies = []

    for entry in conversion_matrix:
      self.ratios.append(entry[0])
      self.efficiencies.append(entry[1])  

  def _convert_direct(self, ratio: float) -> float:
    return ratio / (1 + ratio)

  def _convert_calibrated(self, ratio: float) -> float:
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
  
  def convert(self, ratio : float) -> float:

    if self.calibrated:
      return self._convert_calibrated(ratio)
    else:
      return self._convert_direct(ratio)
    

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