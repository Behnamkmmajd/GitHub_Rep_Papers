"""
Metrics Module
==============

Metric functions for evaluating autoencoder reconstruction quality.
Mirrors the structure of resnet_cae's library_resnet_cae/resnet_cae/metrics.py.
"""

import enum
import numpy as np


class MetricType(enum.Enum):
    MAE = 'MAE'
    MSE = 'MSE'
    LINF = 'LINF'


def compute_mae(input_array, output_array):
    return float(np.mean(np.abs(input_array - output_array)))


def compute_mse(input_array, output_array):
    return float(np.mean((input_array - output_array) ** 2))


def compute_linf(input_array, output_array):
    return float(np.max(np.abs(input_array - output_array)))


_metric_functions = {
    MetricType.MAE: compute_mae,
    MetricType.MSE: compute_mse,
    MetricType.LINF: compute_linf,
}


def compute_metrics_single_array(input_array, output_array, metric_types):
    """Compute a set of metrics for a single pair of arrays.

    Returns:
        dict mapping metric name -> value
    """
    results = {}
    for mt in metric_types:
        results[mt.name] = _metric_functions[mt](input_array, output_array)
    return results
