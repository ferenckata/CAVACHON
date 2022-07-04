
import tensorflow as tf

from cavachon.utils.GeneralUtils import GeneralUtils
from typing import List, Optional, Union

class Prior(tf.keras.layers.Layer):
  """Prior
  Additional parameters used by the CAVACHON model. Including trainable priors pi_y, 
  mean_z_y and logvar_z_y for the latent represnetation in every modalities.

  Attributes:
    n_modalities (int): number of modalities.

    latent_dims (int): number of latent dimensions of each modalities.

    n_components (int): number of components (of VaDE priors) of the priors of 
    latent dimensions of each modalities.

    pi_y (Optional[List[tf.Variable]]): the logits of the marginal probablities for the
    different components of the priors of latent dimensions of each modalities.

    mean_z_y (Optional[List[tf.Variable]]): the means for the different components of 
    the priors of latent dimensions of each modalities.

    logvar_z_y (Optional[List[tf.Variable]]): the log-transformed variancefor the 
    different components of the priors of latent dimensions of each modalities.
  """
  def __init__(
        self,
        prefix: str,
        n_latent_dims: int,
        n_clusters: int):
    super().__init__()
    self.latent_dims: int = n_latent_dims
    self.n_clusters: int = n_clusters   
    self.pi_logit_y = self.add_weight(f'{prefix}:pi_logit_y', (1, n_clusters), trainable=True)
    self.mean_z_y = self.add_weight(f'{prefix}:mean_z_y', (n_latent_dims, n_clusters), trainable=True)
    self.logvar_z_y = self.add_weight(f'{prefix}:logvar_z_y', (n_latent_dims, n_clusters), trainable=True)

    return

  @property
  def pi_y(self):
    return tf.keras.activations.softmax(self.pi_logit_y)

  @property
  def var_z_y(self):
    return tf.math.softplus(self.logvar_z_y)

 