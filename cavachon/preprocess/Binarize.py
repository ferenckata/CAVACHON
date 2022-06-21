from __future__ import annotations
from cavachon.preprocess.PreprocessStep import PreprocessStep

class Binarize(PreprocessStep):

  def __init__(self, name, kwargs):
    super().__init__(name, kwargs)
    
  def execute(self, modality: Modality) -> None:
    threshold = self.kwargs.get('threshold', 1.0)
    modality.adata.X[modality.adata.X >= threshold] = 1.0
    return