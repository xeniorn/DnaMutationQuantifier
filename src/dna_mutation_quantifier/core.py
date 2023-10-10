class Constants:
  """Global constants"""
  peak_search_narrowing_factor = 0.9
  """Sanger peak search range reduction factor, to eliminate crossing over to neighboring peak"""

  default_verbosity = 2
  """for logging"""

  nucleotides = ["A","C","G","T"]
  """Used nucleotide symbols"""

def set_verbosity(verbosity: int):
  """# 0 - mostly silent, 1 - some stuff, 2 - debugging"""
  global global_verbosity
  global_verbosity = verbosity

set_verbosity(Constants.default_verbosity)

def log(message, message_verbosity : int = -1):
  if message_verbosity <= global_verbosity:
    print(message)

