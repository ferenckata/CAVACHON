from __future__ import annotations

import muon as mu
import os
import tensorflow as tf

from cavachon.environment.Constants import Constants
from cavachon.modality.ModalityOrderedMap import ModalityOrderedMap
from cavachon.parser.ConfigParser import ConfigParser
from cavachon.utils.TensorUtils import TensorUtils
from sklearn.preprocessing import LabelEncoder
from typing import Dict, List

class DataLoader:
  """DataLoader
  Data loader to create Tensorflow dataset from MuData.

  Attributes:
    batch_effect_encoder (Dict[str, LabelEncoder]): the encoders used to
    create one-hot encoded batch effect tensor. The keys of the 
    dictionary are formatted as "{modality}:{obs_column}". The 
    LabelEncoder stored the mapping between categorical batch effect
    variables and the numerical representation.

    dataset (tf.data.Dataset): Tensorflow Dataset created from the 
    MuData. Can be used to train/test/validate the CAVACHON model. The 
    field of the dataset includes "{modality}:matrix" (tf.SparseTensor),
    "{modality}:libsize" (tf.Tensor) and "{modality:batch_effect}" 
    (tf.Tensor)

    mdata (mu.MuData): (Single-cell) multi-omics data stored in MuData 
    format.
 
  """
  def __init__(
      self,
      mdata: mu.MuData,
      batch_effect_colnames_dict: Dict[str, List[str]] = dict()) -> None:
    self.batch_effect_encoder: Dict[str, LabelEncoder] = dict()
    self.mdata: mu.MuData = mdata
    self.dataset: tf.data.Dataset = self.create_dataset(batch_effect_colnames_dict)
  
    return

  def create_dataset(
      self,
      batch_effect_colnames_dict: Dict[str, List[str]] = dict()) -> tf.data.Dataset:
    """Create a Tensorflow Dataset based on the MuData provided in the 
    __init__ function.


    Args:
      batch_effect_colnames_dict (Dict[str, List[str]], optional): 
      dictionary of the column name of batch effect in the obs 
      DataFrame, where the keys are the modalities, and values are the 
      batch effect columns (in list) corresponds to the modality. The 
      batch effect columns can be either categorical or continuous. 
      Defaults to None.

    Returns:
      tf.data.Dataset: created Dataset. The field of the dataset 
      includes "{modality}:matrix" (tf.SparseTensor), 
      "{modality}:libsize" (tf.Tensor) and "{modality:batch_effect}" 
      (tf.Tensor)
    """
    field_dict = dict()
    for modality_name in self.mdata.mod.keys():
      adata = self.mdata[modality_name]
      data_tensor = TensorUtils.csr_to_sparse_tensor(adata.X)

      # if batch_effect colname is not specified for the current modality, use zero 
      # matrix as batch effect
      if (batch_effect_colnames_dict is None or modality_name not in batch_effect_colnames_dict):
        batch_effect_tensor = tf.zeros((adata.n_obs, 1))
      else:
        batch_effect_tensor, encoder_dict = TensorUtils.create_tensor_from_df(
            adata.obs, batch_effect_colnames_dict[modality_name]
        )
        for mod, encoder in encoder_dict.items():
          self.batch_effect_encoder[f"{modality_name}:{mod}"] = encoder

      field_dict.setdefault(f"{modality_name}/{Constants.TENSOR_NAME_X}", data_tensor)
      field_dict.setdefault(f"{modality_name}/{Constants.TENSOR_NAME_BATCH}", batch_effect_tensor)
    
    self.dataset = tf.data.Dataset.from_tensor_slices(field_dict)
    return self.dataset

  @classmethod
  def from_modality_ordered_map(cls, modality_map: ModalityOrderedMap) -> DataLoader:
    """Create DataLoader from the dictionary of AnnData.

    Args:
      adata_dict (Dict[str, anndata.AnnData]): dictionary of AnnData, 
      where keys are the modality, values are the corresponding AnnData.

    Returns:
      DataLoader: DataLoader created from the dictionary of AnnData.
    """
    mdata = modality_map.export_mudata()
    mdata.update()
    return cls(mdata)
  
  @classmethod
  def from_h5mu(cls, h5mu_path: str) -> DataLoader:
    """Create DataLoader from h5mu file (of MuData). Note that the 
    different modalities in the MuData needs to be sorted in a way that
    the order of obs DataFrame needs to be the same.

    Args:
      h5mu_path (str): path to the h5mu file.

    Returns:
      DataLoader: DataLoader created from h5mu file.
    """
    path = os.path.realpath(h5mu_path)
    mdata = mu.read(path)
    mdata.update()
    return cls(mdata)

  @classmethod
  def from_config_parser(cls, cp: ConfigParser) -> DataLoader:
    """Create DataLoader from a given config yaml file (config.yaml).

    Args:
      filename (str): the filename of the config yaml.

    Returns:
      DataLoader: DataLoader created from the config file.
    """
    modality_map = ModalityOrderedMap.from_config_parser(cp)
    modality_map.preprocess()
    
    return cls.from_modality_ordered_map(modality_map)
  
  def load_dataset(self, datadir: str) -> None:
    """Load Tensorflow Dataset snapshot.

    Args:
      datadir (str): the data directory of created Tensorflow Dataset 
      snapshot.
    """
    datadir = os.path.realpath(datadir)
    self.dataset = tf.data.experimental.load(datadir)
    return

  def save_dataset(self, datadir: str) -> None:
    """Save Tensorflow Dataset to local storage.

    Args:
      datadir (str): directory where the Tensorflow Dataset snapshot 
      will be save.
    """
    datadir = os.path.realpath(datadir)
    os.makedirs(datadir, exist_ok=True)
    tf.data.experimental.save(self.dataset, datadir)
    return
