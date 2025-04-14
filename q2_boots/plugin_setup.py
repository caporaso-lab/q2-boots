# ----------------------------------------------------------------------------
# Copyright (c) 2024, Caporaso Lab (https://cap-lab.bio).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Range, Collection, Str, Choices, Bool,
                           Float, Metadata, Visualization, Citations, List)

from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence
)
from q2_types.feature_data import (
    FeatureData, Sequence, RNASequence, ProteinSequence
)
from q2_types.sample_data import AlphaDiversity, SampleData

from q2_types.tree import Phylogeny, Rooted

from q2_diversity_lib.alpha import METRICS as alpha_metrics
from q2_diversity_lib.beta import METRICS as beta_metrics
from q2_types.distance_matrix import DistanceMatrix
from q2_types.ordination import PCoAResults

import q2_boots
from q2_boots._examples import (_resample_bootstrap_example,
                                _resample_rarefaction_example,
                                _alpha_rarefaction_example,
                                _alpha_bootstrap_example,
                                _beta_bootstrap_example,
                                _beta_rarefaction_example,
                                _core_metrics_bootstrap_example,
                                _core_metrics_rarefaction_example,
                                _kmer_diversity_bootstrap_example)

citations = Citations.load("citations.bib", package='q2_boots')

plugin = Plugin(
    name='boots',
    version=q2_boots.__version__,
    website='https://github.com/caporaso-lab/q2-boots',
    package='q2_boots',
    short_description=('Bootstrapped and rarefaction-based diversity '
                       'analyses.'),
    description=('A plugin providing bootstrapped and rarefaction-based '
                 '(i.e., resampled) diversity analyses, designed to mirror the '
                 'interface of q2-diversity.'),
    citations=[citations['Raspet2025']]
)


_feature_table_description = 'The input feature table.'
_phylogeny_description = (
    'The phylogenetic tree to use in phylogenetic diversity '
    'calculations. All feature ids in `table` must be present in '
    'this tree, but this tree can contain feature ids that are '
    'not present in `table`.')
_sampling_depth_description = (
    'The total number of observations that each sample in `table` should be '
    'resampled to. Samples where the total number of observations in `table` '
    'is less than `sampling_depth` will be not be included in the output '
    'tables.')
_n_description = 'The number of resampled tables that should be generated.'
_replacement_description = (
    'Resample `table` with replacement (i.e., bootstrap) or without '
    'replacement (i.e., rarefaction).')
_resampled_tables_description = 'The `n` resampled tables.'

# Resampling

_resample_inputs = {
    'table': FeatureTable[Frequency]
}
_resample_parameters = {
    'sampling_depth': Int % Range(1, None),
    'n': Int % Range(1, None),
    'replacement': Bool
}
_resample_outputs = {
    'resampled_tables': Collection[FeatureTable[Frequency]]
}
_resample_input_descriptions = {
    'table': _feature_table_description
}
_resample_parameter_descriptions = {
    'sampling_depth': _sampling_depth_description,
    'n': _n_description,
    'replacement': _replacement_description
}
_resample_output_descriptions = {
    'resampled_tables': _resampled_tables_description
}

plugin.pipelines.register_function(
    function=q2_boots.resample,
    inputs=_resample_inputs,
    parameters=_resample_parameters,
    outputs=_resample_outputs,
    input_descriptions=_resample_input_descriptions,
    parameter_descriptions=_resample_parameter_descriptions,
    output_descriptions=_resample_output_descriptions,
    name='Resample feature table, returning `n` feature tables.',
    description=('Resample `table` to `sampling_depth` total observations with '
                 'replacement (i.e., bootstrapping) or without replacement '
                 '(i.e., rarefaction) `n` times, to generate `n` resampled '
                 'feature tables.'),
    examples={
        'Generate 10 bootstrapped tables.': _resample_bootstrap_example,
        'Generate 10 rarefied tables.': _resample_rarefaction_example
    }
)

_diversity_inputs = {
    'table': FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
    'phylogeny': Phylogeny[Rooted]
}

_diversity_input_descriptions = {
    'table': _feature_table_description,
    'phylogeny': _phylogeny_description
}

_alpha_average_parameters = {
    'average_method': Str % Choices('mean', 'median')
}

_alpha_average_parameter_descriptions = {
    'average_method': 'Method to use for averaging.'
}

_average_alpha_diversity_description = (
    'The average alpha diversity vector.')

plugin.methods.register_function(
    function=q2_boots.alpha_average,
    inputs={
        'data': Collection[SampleData[AlphaDiversity]]
    },
    parameters=_alpha_average_parameters,
    outputs={
        'average_alpha_diversity': SampleData[AlphaDiversity]
    },
    input_descriptions={
        'data': 'Alpha diversity vectors to be averaged.'
    },
    output_descriptions={
        'average_alpha_diversity': _average_alpha_diversity_description

    },
    parameter_descriptions=_alpha_average_parameter_descriptions,
    name='Average alpha diversity vectors.',
    description=('Compute the per-sample average across a collection of alpha '
                 'diversity vectors computed from the same samples.')
)

_alpha_collection_parameters = {
    'sampling_depth': Int % Range(1, None),
    'metric': Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                            alpha_metrics['NONPHYLO']['UNIMPL'] |
                            alpha_metrics['PHYLO']['IMPL'] |
                            alpha_metrics['PHYLO']['UNIMPL']),
    'n': Int % Range(1, None),
    'replacement': Bool
}

_alpha_collection_parameter_descriptions = {
    'sampling_depth': _sampling_depth_description,
    'metric': 'The alpha diversity metric to be computed.',
    'n': _n_description,
    'replacement': _replacement_description
}

plugin.pipelines.register_function(
    function=q2_boots.alpha_collection,
    inputs=_diversity_inputs,
    parameters=_alpha_collection_parameters,
    outputs={'alpha_diversities': Collection[SampleData[AlphaDiversity]]},
    input_descriptions=_diversity_input_descriptions,
    parameter_descriptions=_alpha_collection_parameter_descriptions,
    output_descriptions={
        'alpha_diversities': ('`n` alpha diversity vectors, each containing '
                              'per-sample alpha diversity scores for the same '
                              'samples.'),
    },
    name='Perform resampled alpha diversity, returning `n` result vectors.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes the specified '
                 'alpha diversity metric on each resulting `table`. The '
                 'resulting artifacts can be used, for example, to explore the '
                 'variance across `n` iterations of resampling.')
)

_alpha_parameters = _alpha_collection_parameters | _alpha_average_parameters
_alpha_parameter_descriptions = (_alpha_collection_parameter_descriptions |
                                 _alpha_average_parameter_descriptions)

plugin.pipelines.register_function(
    function=q2_boots.alpha,
    inputs=_diversity_inputs,
    parameters=_alpha_parameters,
    outputs={'average_alpha_diversity': SampleData[AlphaDiversity]},
    input_descriptions=_diversity_input_descriptions,
    parameter_descriptions=_alpha_parameter_descriptions,
    output_descriptions={
        'average_alpha_diversity': _average_alpha_diversity_description,
    },
    name='Perform resampled alpha diversity, returning average result vector.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes the specified '
                 'alpha diversity metric on each resulting `table`. The '
                 'resulting artifacts are then averaged using the method '
                 'specified by `average_method`, and the resulting average '
                 'per-sample alpha diversities are returned.'),
    examples={
        'Bootstrapped observed features.': _alpha_bootstrap_example,
        'Rarefaction-based observed features.': _alpha_rarefaction_example
    }
)

_beta_average_parameters = {
    'average_method': Str % Choices(['non-metric-mean',
                                     'non-metric-median',
                                     'medoid'])
}

_beta_average_parameter_descriptions = {
    'average_method': 'Method to use for averaging.'
}

plugin.methods.register_function(
    function=q2_boots.beta_average,
    inputs={
        'data': Collection[DistanceMatrix],
    },
    parameters=_beta_average_parameters,
    outputs={'average_distance_matrix': DistanceMatrix},
    input_descriptions={
        'data': 'Distance matrices to be average.'
    },
    output_descriptions={
        'average_distance_matrix': 'The average distance matrix.',
    },
    parameter_descriptions=_beta_average_parameter_descriptions,
    name='Average beta diversity distance matrices.',
    description=('Compute the average distance matrix across a collection of '
                 'distance matrices.')
)

_beta_collection_parameters = {
                'metric': Str % Choices(beta_metrics['NONPHYLO']['IMPL'] |
                                        beta_metrics['NONPHYLO']['UNIMPL'] |
                                        beta_metrics['PHYLO']['IMPL'] |
                                        beta_metrics['PHYLO']['UNIMPL']),
                'pseudocount': Int % Range(1, None),
                'replacement': Bool,
                'n': Int % Range(1, None),
                'sampling_depth': Int % Range(1, None),
                'bypass_tips': Bool,
                'variance_adjusted': Bool,
                'alpha': Float % Range(0, 1, inclusive_end=True)
}

_beta_collection_parameter_descriptions = {
    'metric': 'The beta diversity metric to be computed.',
    'pseudocount': ('A pseudocount to handle zeros for compositional '
                    'metrics.  This is ignored for other metrics.'),
    'replacement': _replacement_description,
    'n': _n_description,
    'sampling_depth': _sampling_depth_description,
    'bypass_tips': ('Ignore tips of tree in phylogenetic diversity '
                    'calculations, trading specificity for reduced compute '
                    'time. This has the effect of collapsing the phylogeny, '
                    'and is analogous (in concept) to moving from 99% to 97% '
                    'OTUs.'),
    'variance_adjusted': ('Perform variance adjustment based on Chang et al. '
                          'BMC Bioinformatics (2011) for phylogenetic '
                          'diversity metrics.'),
    'alpha': ('The alpha value used with the generalized UniFrac metric.')
}


plugin.pipelines.register_function(
    function=q2_boots.beta_collection,
    inputs=_diversity_inputs,
    parameters=_beta_collection_parameters,
    outputs={'distance_matrices': Collection[DistanceMatrix]},
    input_descriptions=_diversity_input_descriptions,
    output_descriptions={
        'distance_matrices': ('`n` beta diversity distance matrices, each '
                              'containing distances between all pairs of '
                              'samples and computed from resampled feature '
                              'tables.')
    },
    parameter_descriptions=_beta_collection_parameter_descriptions,
    name='Perform resampled beta diversity, returning `n` distance matrices.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes the specified '
                 'beta diversity metric on each resulting `table`. The '
                 'resulting artifacts can be used, for example, to explore the '
                 'variance across `n` iterations of resampling.'),
    examples={
        'Bootstrapped Bray-Curtis.': _beta_bootstrap_example,
        'Rarefaction-based Bray-Curtis.': _beta_rarefaction_example
    }
)

_beta_parameters = _beta_collection_parameters | _beta_average_parameters
_beta_parameter_descriptions = (_beta_collection_parameter_descriptions |
                                _beta_average_parameter_descriptions)

plugin.pipelines.register_function(
    function=q2_boots.beta,
    inputs=_diversity_inputs,
    parameters=_beta_parameters,
    outputs=[('average_distance_matrix', DistanceMatrix)],
    input_descriptions=_diversity_input_descriptions,
    parameter_descriptions=_beta_parameter_descriptions,
    output_descriptions={
        'average_distance_matrix': 'The average distance matrix.'},
    name='Perform resampled beta diversity, returning average distance matrix.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes the specified '
                 'beta diversity metric on each resulting `table`. The '
                 'resulting artifacts are then averaged using the method '
                 'specified by `average_method`, and the resulting average '
                 'beta diversity distance matrix is returned.')
)

plugin.pipelines.register_function(
    function=q2_boots.core_metrics,
    inputs=_diversity_inputs,
    parameters={
        'metadata': Metadata,
        'n': Int % Range(1, None),
        'sampling_depth': Int % Range(1, None),
        'alpha_average_method': Str % Choices('mean', 'median'),
        'beta_average_method': Str % Choices('non-metric-mean',
                                             'non-metric-median',
                                             'medoid'),
        'replacement': Bool
    },
    outputs=[
        ('resampled_tables', Collection[FeatureTable[Frequency]]),
        ('alpha_diversities', Collection[SampleData[AlphaDiversity]]),
        ('distance_matrices', Collection[DistanceMatrix]),
        ('pcoas', Collection[PCoAResults]),
        ('emperor_plots', Collection[Visualization]),
    ],
    input_descriptions=_diversity_input_descriptions,
    parameter_descriptions={
        'metadata': 'The sample metadata used in generating Emperor plots.',
        'n': _n_description,
        'sampling_depth': _sampling_depth_description,
        'alpha_average_method': 'Method to use for averaging alpha diversity.',
        'beta_average_method': 'Method to use for averaging beta diversity.',
        'replacement': _replacement_description
    },
    output_descriptions={
        'resampled_tables': _resampled_tables_description,
        'alpha_diversities': 'Average alpha diversity vector for each metric.',
        'distance_matrices': ('Average beta diversity distance matrix for '
                              'each metric.'),
        'pcoas': 'PCoA matrix for each beta diversity metric.',
        'emperor_plots': 'Emperor plot for each beta diversity metric.'
    },
    name='Perform resampled "core metrics" analysis.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample, and then computes common alpha '
                 'and beta diversity on each resulting `table`. The '
                 'resulting artifacts are then averaged using the method '
                 'specified by `alpha_average_method` and '
                 '`beta_average_method` parameters. The resulting average '
                 'alpha and beta diversity artifacts are returned, along with '
                 'PCoA matrices and Emperor plots.'),
    examples={
        'Bootstrapped core metrics.': _core_metrics_bootstrap_example,
        'Rarefaction-based core metrics.': _core_metrics_rarefaction_example
    }
)

plugin.pipelines.register_function(
    function=q2_boots.kmer_diversity,
    inputs={'table': FeatureTable[Frequency |
                                  RelativeFrequency |
                                  PresenceAbsence],
            'sequences': FeatureData[Sequence |
                                     RNASequence |
                                     ProteinSequence]},
    parameters={
        'metadata': Metadata,
        'n': Int % Range(1, None),
        'sampling_depth': Int % Range(1, None),
        'alpha_metrics':
            List[Str % Choices(alpha_metrics['NONPHYLO']['IMPL'] |
                               alpha_metrics['NONPHYLO']['UNIMPL'])],
        'beta_metrics': List[Str % Choices(
                                beta_metrics['NONPHYLO']['IMPL'] |
                                beta_metrics['NONPHYLO']['UNIMPL'])],
        'alpha_average_method': Str % Choices('mean', 'median'),
        'beta_average_method': Str % Choices('non-metric-mean',
                                             'non-metric-median',
                                             'medoid'),
        'replacement': Bool,
        'kmer_size': Int,
        'tfidf': Bool,
        'max_df': Float % Range(0, 1, inclusive_start=True,
                                inclusive_end=True) | Int,
        'min_df': Float % Range(0, 1, inclusive_start=True,
                                inclusive_end=False) | Int,
        'max_features': Int,
        'norm': Str % Choices(['None', 'l1', 'l2']),
        'pc_dimensions': Int,
        'color_by': Str
    },
    outputs=[
        ('resampled_tables', Collection[FeatureTable[Frequency]]),
        ('kmer_tables', Collection[FeatureTable[Frequency]]),
        ('alpha_diversities', Collection[SampleData[AlphaDiversity]]),
        ('distance_matrices', Collection[DistanceMatrix]),
        ('pcoas', Collection[PCoAResults]),
        ('scatter_plot', Visualization),
    ],
    input_descriptions={'table': _feature_table_description,
                        'sequences': "Input sequences for kmerization."},
    parameter_descriptions={
        'metadata': 'The sample metadata used in generating Emperor plots.',
        'n': _n_description,
        'sampling_depth': _sampling_depth_description,
        'alpha_average_method': 'Method to use for averaging alpha diversity.',
        'beta_average_method': 'Method to use for averaging beta diversity.',
        'replacement': _replacement_description,
        'kmer_size': 'Length of kmers to generate.',
        'tfidf': 'If True, kmers will be scored using TF-IDF and output '
                 'frequencies will be weighted by scores. If False, kmers are '
                 'counted without TF-IDF scores.',
        'max_df': 'Ignore kmers that have a frequency strictly higher than '
                  'the given threshold. If float, the parameter represents a '
                  'proportion of sequences, if an integer it represents an '
                  'absolute count.',
        'min_df': 'Ignore kmers that have a frequency strictly lower than '
                  'the given threshold. If float, the parameter represents a '
                  'proportion of sequences, if an integer it represents an '
                  'absolute count.',
        'max_features': 'If not None, build a vocabulary that only considers '
                        'the top max_features ordered by frequency (or TF-IDF '
                        'score).',
        'norm': 'Normalization procedure applied to TF-IDF scores. Ignored '
                'if tfidf=False. l2: Sum of squares of vector elements is 1. '
                'l1: Sum of absolute values of vector elements is 1.',
        'pc_dimensions': 'Number of principal coordinate dimensions to keep '
                         'for plotting.',
        'color_by': 'Categorical measure from the input Metadata that '
                    'should be used for color-coding the scatterplot.'
    },
    output_descriptions={
        'resampled_tables': _resampled_tables_description,
        'kmer_tables': 'The `n` kmer tables.',
        'alpha_diversities': 'Average alpha diversity vector for each metric.',
        'distance_matrices': ('Average beta diversity distance matrix for '
                              'each metric.'),
        'pcoas': 'PCoA matrix for each beta diversity metric.',
        'scatter_plot': ('Scatter plot including alpha diversity and '
                         'pcoa results for all selected metrics.')
    },
    name='Perform resampled "core metrics" analysis on kmerized features.',
    description=('Given a single feature table as input, this action resamples '
                 'the feature table `n` times to a total frequency of '
                 '`sampling depth` per sample. It then splits all input '
                 'sequences into kmers, and computes common alpha '
                 'and beta diversity on each resulting kmer table. The '
                 'resulting artifacts are then averaged using the method '
                 'specified by `alpha_average_method` and '
                 '`beta_average_method` parameters. The resulting average '
                 'alpha and beta diversity artifacts are returned, along with '
                 'a scatter plot integrated all alpha diversity metrics and '
                 'the PCoA axes for all beta diversity metrics.'),
    examples={'Bootstrapped kmer diversity': _kmer_diversity_bootstrap_example},
    citations=[citations['Bokulich2024']]
)
