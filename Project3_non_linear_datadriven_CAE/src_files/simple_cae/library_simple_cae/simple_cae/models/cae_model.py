"""
Convolutional Autoencoder Model
================================

Architecture (from Wang et al. 2024):
    Encoder: 256 -> 128 -> 64 -> 32 -> 16 -> 8 -> 4 -> latent
    Decoder: latent -> 4 -> 8 -> 16 -> 32 -> 64 -> 128 -> 256
"""

import torch
import torch.nn as nn


class ConvolutionalAutoencoder(nn.Module):
    """
    Convolutional Autoencoder for compressing 256x256 snapshots to latent space.

    Architecture:
        Encoder: Conv2d+ReLU+MaxPool2d x6 -> Flatten -> Linear
        Decoder: Linear -> Reshape -> (Upsample+Conv2d+ReLU) x6

    Parameters
    ----------
    latent_dim : int
        Size of the latent vector.
    in_channels : int
        Number of input channels (1 for single field, N for multi-channel).
    """

    def __init__(self, latent_dim=4, in_channels=1):
        super(ConvolutionalAutoencoder, self).__init__()
        self.latent_dim = latent_dim
        self.in_channels = in_channels

        # ENCODER: 256x256 x in_channels -> latent_dim
        self.encoder = nn.Sequential(
            nn.Conv2d(in_channels, 16, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),

            nn.Conv2d(16, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),

            nn.Conv2d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),

            nn.Conv2d(64, 128, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),

            nn.Conv2d(128, 256, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),

            nn.Conv2d(256, 512, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2),
        )

        self.fc_encode = nn.Linear(4 * 4 * 512, latent_dim)

        # DECODER: latent_dim -> 256x256 x in_channels
        self.fc_decode = nn.Linear(latent_dim, 4 * 4 * 512)

        self.decoder = nn.Sequential(
            nn.Upsample(scale_factor=2),
            nn.Conv2d(512, 256, kernel_size=3, padding=1),
            nn.ReLU(),

            nn.Upsample(scale_factor=2),
            nn.Conv2d(256, 128, kernel_size=3, padding=1),
            nn.ReLU(),

            nn.Upsample(scale_factor=2),
            nn.Conv2d(128, 64, kernel_size=3, padding=1),
            nn.ReLU(),

            nn.Upsample(scale_factor=2),
            nn.Conv2d(64, 32, kernel_size=3, padding=1),
            nn.ReLU(),

            nn.Upsample(scale_factor=2),
            nn.Conv2d(32, 16, kernel_size=3, padding=1),
            nn.ReLU(),

            nn.Upsample(scale_factor=2),
            nn.Conv2d(16, in_channels, kernel_size=3, padding=1),
        )

    def encode(self, x):
        x = self.encoder(x)
        x = x.view(x.size(0), -1)
        z = self.fc_encode(x)
        return z

    def decode(self, z):
        x = self.fc_decode(z)
        x = x.view(x.size(0), 512, 4, 4)
        x = self.decoder(x)
        return x

    def forward(self, x):
        z = self.encode(x)
        x_reconstructed = self.decode(z)
        return x_reconstructed
