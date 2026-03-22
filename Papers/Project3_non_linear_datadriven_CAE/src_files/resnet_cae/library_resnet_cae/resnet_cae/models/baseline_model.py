import logging
from functools import partial

import torch.nn as nn

from resnet_cae.models.blocks import Upsample, Downsample, BasicBlock, ResnetBlock

logger = logging.getLogger(__name__)


class Encoder(nn.Module):
    def __init__(self,
                 dim: int,
                 dim_mults: tuple,
                 channels: int,
                 z_channels: int,
                 block_type: int,
                 resnet_block_groups: int = 4,
                 act_type: str = 'silu'
                 ):
        super().__init__()
        self.channels = channels
        self.init_conv = nn.Conv2d(channels, dim, 1, padding=0)

        dims = list(map(lambda m: dim * m, dim_mults))
        in_out = list(zip(dims[:-1], dims[1:]))

        if block_type == 0:
            block_class = BasicBlock
        elif block_type == 1:
            block_class = partial(ResnetBlock, groups=resnet_block_groups, act_type=act_type)
        else:
            raise ValueError(f"Invalid block type: {block_type}")

        self.downs = nn.ModuleList([])

        for ind, (dim_in, dim_out) in enumerate(in_out):
            is_last = ind == len(in_out) - 1

            self.downs.append(block_class(dim_in, dim_in))
            self.downs.append(block_class(dim_in, dim_in))
            self.downs.append(Downsample(dim_in, dim_out) if not is_last else nn.Conv2d(dim_in, dim_out, 3, padding=1))

        mid_dim = dims[-1]
        self.mid_block_1 = block_class(mid_dim, mid_dim)
        self.mid_block_2 = block_class(mid_dim, mid_dim)
        self.final_block = nn.Conv2d(mid_dim, z_channels, kernel_size=3, stride=1, padding=1)

    def forward(self, x):
        h = self.init_conv(x)

        for block in self.downs:
            h = block(h)

        h = self.mid_block_1(h)
        h = self.mid_block_2(h)
        h = self.final_block(h)

        if not hasattr(self, "_printed_latent"):
            print("LATENT SPACE SHAPE =", h.shape)
            self._printed_latent = True  # ensures it prints only once

        return h


class Decoder(nn.Module):
    """Reconstruct input x from latent variable z
    """

    def __init__(self,
                 dim: int,
                 dim_mults: tuple,
                 channels: int,
                 z_channels: int,
                 block_type: int,
                 resnet_block_groups: int = 4,
                 final_kernel_size: int = 3,
                 act_type: str = 'silu'
                 ):
        super().__init__()
        self.channels = channels

        dims = list(map(lambda m: dim * m, dim_mults))
        in_out = list(zip(dims[:-1], dims[1:]))

        if block_type == 0:
            block_class = BasicBlock
        elif block_type == 1:
            block_class = partial(ResnetBlock, groups=resnet_block_groups, act_type=act_type)
        else:
            raise ValueError(f"Invalid block type: {block_type}")

        self.ups = nn.ModuleList([])

        self.init_conv = nn.Conv2d(z_channels, dims[-1], 3, padding=1)

        for ind, (dim_in, dim_out) in enumerate(reversed(in_out)):
            is_last = ind == len(in_out) - 1

            self.ups.append(block_class(dim_out, dim_out))
            self.ups.append(block_class(dim_out, dim_out))
            self.ups.append(Upsample(dim_out, dim_in) if not is_last else nn.Conv2d(dim_out, dim_in, 3, padding=1))

        out_dim = dims[0]
        self.mid_block_1 = block_class(out_dim, out_dim)
        self.mid_block_2 = block_class(out_dim, out_dim)

        if final_kernel_size == 3:
            self.final_block = nn.Conv2d(out_dim, channels, kernel_size=3, stride=1, padding=1)
        elif final_kernel_size == 1:
            self.final_block = nn.Conv2d(out_dim, channels, kernel_size=1, stride=1, padding=0)
        else:
            raise ValueError(f"Invalid final kernel size: {final_kernel_size}")

    def forward(self, z):
        h = self.init_conv(z)

        for block in self.ups:
            h = block(h)

        h = self.mid_block_1(h)
        h = self.mid_block_2(h)
        h = self.final_block(h)
        return h


class ConvAutoencoderBaseline(nn.Module):
    def __init__(self,
                 dim,
                 dim_mults,
                 channels,
                 z_channels,
                 block_type,
                 act_type='silu',
                 ):
        super().__init__()
        self.encoder = Encoder(dim, dim_mults, channels, z_channels, block_type=block_type, act_type=act_type)
        self.decoder = Decoder(dim, dim_mults, channels, z_channels, block_type=block_type, act_type=act_type)

        num_params = sum(p.numel() for p in self.parameters())
        logger.info(f'Constructed ConvAutoencoderBaseline with {num_params} parameters [act_type = {act_type}]')
        logger.debug(f'Model architecture: \n{self}')

    def encode(self, x):
        z = self.encoder(x)
        return z

    def decode(self, z):
        dec = self.decoder(z)
        return dec

    def forward(self, input):
        z = self.encode(input)
        x = self.decode(z)
        return x
