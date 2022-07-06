from __future__ import annotations
from cavachon.filter.FilterStep import FilterStep
import scanpy

class FilterCells(FilterStep):

  def __init__(self, name, args):
    super().__init__(name, args)

  def execute(self, modality: Modality) -> None:
    self.kwargs['inplace'] = True
    scanpy.pp.filter_cells(modality.adata, **self.kwargs)
    modality.n_obs = modality.adata.n_obs
    modality.n_vars = modality.adata.n_vars
    return