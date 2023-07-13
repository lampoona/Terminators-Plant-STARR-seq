import torch
import torch.nn as nn
from eugene.models.base import BaseModel
from typing import Union, Callable
from eugene.models.base import BasicConv1D, BasicFullyConnectedModule
from pytorch_lightning.utilities.cli import LightningArgumentParser

# Define a DenseNet model:
# adapted from https://doi.org/10.1016/j.xplc.2022.100455
# pytorch implementation inspired by https://uvadlc-notebooks.readthedocs.io/en/latest/tutorial_notebooks/tutorial5/Inception_ResNet_DenseNet.html

# Define a registry of activation functions (copied from the github version of EUGENe)
ACTIVATION_REGISTRY = {
    "relu": nn.ReLU,
    "leaky_relu": nn.LeakyReLU,
    "gelu": nn.GELU,
    "elu": nn.ELU,
    "sigmoid": nn.Sigmoid,
    "tanh": nn.Tanh,
    "softplus": nn.Softplus
}

# Define a dense layer (applies two succesive convolutional layers to the input and concatenates the results with the input):
class DenseLayer(nn.Module):
    """
    Inputs:
        c_in - number of input channels
        growth_rate - number of final ouput channels
        R - multiplier of growth_rate to determine the number of ouput channels in the first convolutional layer
        kernel_size = size of the convolving kernel
        act_fn = activation function used for each convolutional layer
        drop_rate = probability for the dropout layer
    """
    def __init__(
        self,
        c_in,
        growth_rate = 4,
        R = 1,
        kernel_size = 3,
        act_fn = nn.ReLU,
        drop_rate = 0.2
    ):
        super().__init__()
        net = [
            nn.Conv1d(in_channels = c_in, out_channels = growth_rate * R, kernel_size = kernel_size, stride = 1, padding = 'same'),
            act_fn(),
            nn.Conv1d(in_channels = growth_rate * R, out_channels = growth_rate, kernel_size = kernel_size, stride = 1, padding = 'same'),
            act_fn()
        ]
        if drop_rate != 0.0:
            net.append(nn.Dropout(p = drop_rate))
        self.net = nn.Sequential(*net)

    def forward(self, x):
        out = self.net(x)
        out = torch.cat([x, out], dim = 1)
        return out

# Define a dense block (combines multiple dense layers):
class DenseBlock(nn.Module):
    """
    Inputs:
        c_in - number of input channels
        num_layer - number of dense layers
        growth_rate - number of final ouput channels
        R - multiplier of growth_rate to determine the number of ouput channels in the first convolutional layer
        kerner_size = size of the convolving kernel
        act_fn = activation function used for each convolutional layer
        drop_rate = probability for the dropout layer
    """
    def __init__(
        self,
        c_in,
        num_layers = 4,
        growth_rate = 4,
        R = 1,
        kernel_size = 3,
        act_fn = nn.ReLU,
        drop_rate = 0.2
    ):
        super().__init__()
        layers = []
        for layer_idx in range(num_layers):
            layers.append(
                DenseLayer(
                    c_in = c_in + layer_idx * growth_rate, # input channels are original plus the feature maps from previous layers
                    growth_rate = growth_rate,
                    R = R,
                    kernel_size = kernel_size,
                    act_fn = act_fn,
                    drop_rate = drop_rate
                )
            )
        self.block = nn.Sequential(*layers)

    def forward(self, x):
        out = self.block(x)
        return out

# Define a transition layer (compresses the output from the previous dense block):
class TransitionLayer(nn.Module):
    """
    Inputs:
        c_in - number of input channels,
        compression_factor = factor by which the output channels from the previous dense block are compressed
        pool_size = kernel size for the average pooling layer
        act_fn = activation function used for the convolutional layer
    """
    def __init__(
        self,
        c_in,
        compression_factor = 2,
        pool_size = 2,
        act_fn = nn.ReLU
    ):
        super().__init__()
        self.transition = nn.Sequential(
            nn.Conv1d(in_channels = c_in, out_channels = c_in // compression_factor, kernel_size = 1, stride = 1, padding = 'same'),
            act_fn(),
            nn.AvgPool1d(kernel_size = pool_size)
        )

    def forward(self, x):
        return self.transition(x)

# Define a DenseNet module:
class BasicDenseNet(nn.Module):
    """Generates a PyTorch module with a basic DenseNet architecture
    (i.e. groups of dense layers separated by transition layers)
    Parameters
    ----------
    input_dim : tuple
        input dimensions
    layers : list
        number of layers for each dense block
    growth_rates : list-like or int
        growth rate for each dense block (i.e. number of channels added in each dense layer)
    R : list-like or int
        multiplier of growth_rate to determine the number of ouput channels in the first convolutional layer of each dense layer
    conv_kernels : list-like or int
        conv kernel size for each conv layer
    activation : pytorch activation function or str
        activation function used after all convolutional layers
    dropout_rates : list-like or float
        dropout rates for each convolutional layer
    compression : list-like or int
        factor by which each transition layers compresses the channels from the previous layer
    pool_kernels : list-like or int
        average pooling kernel size for each transition layer
    omit_final_transition : boolean
        omit a transition layer after the last dense block
    Returns
    -------
    nn.Module
        The instantiated DenseNet module.
    """
    def __init__(
        self,
        input_dim: tuple,
        layers: list,
        growth_rates: Union[list, int],
        R: Union[list, int],
        conv_kernels: Union[list, int] = 3,
        dropout_rates: Union[list, float] = 0.2,
        compression: Union[list, int] = 2,
        pool_kernels: Union[list, int] = 2,
        activation: Union[str, Callable] = "relu",
        omit_final_transition: bool = True
    ):
        super(BasicDenseNet, self).__init__()
        assert isinstance(layers, list)
        growth_rates = growth_rates if isinstance(growth_rates, list) else [growth_rates] * len(layers)
        R = R if isinstance(R, list) else [R] * len(layers)
        conv_kernels = conv_kernels if isinstance(conv_kernels, list) else [conv_kernels] * len(layers)
        dropout_rates = dropout_rates if isinstance(dropout_rates, list) else [dropout_rates] * len(layers)
        compression = compression if isinstance(compression, list) else [compression] * len(layers)
        pool_kernels = pool_kernels if isinstance(pool_kernels, list) else [pool_kernels] * len(layers)
        activation = ACTIVATION_REGISTRY[activation]

        net = []
        c_in, len_in = input_dim
        for i in range(len(layers)):
            net.append(
                DenseBlock(
                    c_in = c_in,
                    num_layers = layers[i],
                    growth_rate = growth_rates[i],
                    R = R[i],
                    kernel_size = conv_kernels[i],
                    act_fn = activation,
                    drop_rate = dropout_rates[i]
                )
            )
            c_in = c_in + layers[i] * growth_rates[i]
            if i < (len(layers) - 1) or not omit_final_transition:
                net.append(
                    TransitionLayer(
                        c_in = c_in,
                        compression_factor = compression[i],
                        pool_size = pool_kernels[i],
                        act_fn = activation
                    )
                )
                c_in = c_in // compression[i]
                len_in = len_in // pool_kernels[i]
        
        self.module = nn.Sequential(*net)
        self.out_channels = c_in
        self.flatten_dim = int(c_in * len_in)

    def forward(self, x):
        return self.module(x)

# Define a function to build the DenseNet model:
class DenseNet(BaseModel):
    def __init__(
        self,
        input_len: int,
        output_dim: int,
        conv_kwargs: dict,
        strand: str = "ss",
        task: str = "regression",
        aggr: str = None,
        loss_fxn: str = "mse",
        dn_kwargs: dict = {},
        fc_kwargs: dict = {},
        **kwargs
    ):
        super().__init__(
            input_len, 
            output_dim, 
            strand = strand, 
            task = task, 
            aggr = aggr, 
            loss_fxn = loss_fxn,
            **kwargs
        )
        self.convnet = BasicConv1D(
            input_len = input_len, 
            **conv_kwargs
        )
        self.densenet = BasicDenseNet(
            input_dim = (self.convnet.out_channels, 170), # self.convnet.flatten_dim doesn't work with padding = 'same' in this version; fixed in the github version; use the line below with the fixed code
            # input_dim = (self.convnet.out_channels, self.convnet.flatten_dim / self.convnet.out_channels),
            **dn_kwargs
        )
        self.fcnet = BasicFullyConnectedModule(
            input_dim = self.densenet.flatten_dim,
            output_dim = output_dim, 
            **fc_kwargs
        )

    def forward(self, x, x_rev_comp = None):
        x = self.convnet(x)
        x = self.densenet(x)
        x = x.view(x.size(0), self.densenet.flatten_dim)
        x = self.fcnet(x)
        return x

    def load_from_config(config_file):
        parser = LightningArgumentParser()
        model_type = DenseNet
        parser.add_lightning_class_args(model_type, nested_key = 'model')
        model_yml = parser.parse_path(config_file)
        return model_type(**model_yml['model'])
