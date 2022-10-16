from cavachon.config.OptimizerConfig import OptimizerConfig
from cavachon.config.ConfigMapping import ConfigMapping
from typing import Any, Mapping

class TrainingConfig(ConfigMapping):
  """TrainingConfig

  Config for training.

  Attributes
  ----------
  optimizer: OptimizerConfig
      config for optimizer.
  
  max_n_epochs: int
      maximum number of epochs for training.
      
  """
  def __init__(self, config: Mapping[str, Any]):
    """Constructor for TrainingConfig. 

    Parameters
    ----------
    config: Mapping[str, Any]:
        training config in mapping format.
    
    """
    self.optimizer: OptimizerConfig
    self.max_n_epochs: int
    super().__init__(config)