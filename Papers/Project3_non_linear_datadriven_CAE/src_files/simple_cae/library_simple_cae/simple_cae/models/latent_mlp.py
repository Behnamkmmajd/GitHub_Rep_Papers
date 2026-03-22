"""
Latent MLP Model
================

Multi-Layer Perceptron to map parameters (mu, t) -> latent codes z.
This is the second network in CNN-ROM (after the CAE).

Architecture:
    Linear(2, 64) -> ReLU -> Linear(64, 64) -> ReLU -> Linear(64, latent_dim)
"""

import torch
import torch.nn as nn


class LatentMLP(nn.Module):
    """
    MLP to map normalized (mu, t) -> latent codes z.
    """

    def __init__(self, input_dim=2, hidden_dim=64, latent_dim=4):
        super(LatentMLP, self).__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim

        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, latent_dim),
        )

    def forward(self, x):
        return self.network(x)
