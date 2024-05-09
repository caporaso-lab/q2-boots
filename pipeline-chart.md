```mermaid
flowchart TD
    
    feature-table[FeatureTable] --> resample
    n(["n (Int)"]) --> resample
    sampling-depth(["sampling-depth (Int)"]) --> resample
    phylogeny-alpha[phylogeny] -- optional --> div-alpha
    alpha-metric(["metric (alpha, String)"]) --> div-alpha
    resample{qiime feature-table resample} -->
    feature-collection[[FeatureTable]] -- looped call --> div-alpha
    feature-collection -- looped call --> beta-diversity
    div-alpha{qiime diversity alpha} --> sd-collection
    alpha-rep(["representative (alpha, String)"]) --> alpha-avg
    sd-collection[["SampleData[AlphaDiversity]"]] -->
    alpha-avg{qiime boots alpha-average} -->
    I["SampleData[AlphaDiversity]"]

    beta-diversity{qiime diversity beta} --> dms
    beta-metric(["metric (beta, String)"]) --> beta-diversity
    phylogeny-beta[phylogeny] -- optional --> beta-diversity

    dms[[DistanceMatrix]] --> beta-average
    beta-rep(["representative (beta, String)"]) -->
    beta-average{qiime boots beta-average} -->
    dm[DistanceMatrix]
```
  