## Introduction

This project investigates the relationship between personality characteristics and patterns of psychoactive-drug consumption using statistical-learning methods.

The dataset contains survey responses from 1,885 individuals collected between March 2011 and March 2012. It includes demographic characteristics, personality measurements, impulsivity and sensation-seeking scores, and self-reported consumption frequencies for 18 central-nervous-system psychoactive substances.

The analysis addresses three related questions:

- Can meaningful personality profiles be identified from psychological measurements?
- Can multidimensional drug-consumption patterns be transformed into an interpretable response variable for classification?
- Can drug-consumption profiles be predicted from personality characteristics?

### Personality profiling

Personality profiles are first identified from neuroticism, extraversion, openness, agreeableness, conscientiousness, impulsivity, and sensation-seeking measurements.

Hierarchical clustering is used to explore the number and structure of personality groups. Based on this analysis, three clusters are selected, and their centroids are subsequently used to initialize k-means clustering. The resulting profiles represent distinct combinations of personality characteristics, ranging from more cautious and disciplined individuals to sensation-seeking and emotionally unstable profiles.

### Drug-consumption profiling

Drug consumption is represented by multiple correlated variables describing the recency or frequency of use of different substances. Principal Component Analysis (PCA) is therefore used to reduce the dimensionality of the drug-consumption space.

Several criteria are considered when determining the number of components to retain, including:

- explained variance and scree analysis;
- the Kaiser criterion;
- Horn's parallel analysis;
- Velicer's MAP test;
- Very Simple Structure (VSS).

Seven components are retained and subjected to Varimax rotation to improve interpretation. Together, these rotated components explain approximately 67% of the total variation in the drug-consumption variables.

To reflect differences in the dependence potential of different substances, external physical- and psychological-dependence scores are incorporated into the analysis. These scores are propagated through the component loadings to construct dependency-based weights for the rotated component space.

A weighted Manhattan distance is then defined and used for clustering individuals according to their drug-consumption patterns. Partitioning Around Medoids (PAM) and hierarchical clustering are compared using elbow analysis, silhouette widths, cluster sizes, and interpretability.

Although some hierarchical solutions achieve higher average silhouette widths, they produce highly imbalanced clusters. PAM with three clusters provides a more useful partition for identifying distinct drug-consumption profiles and is therefore used to define the response variable for the subsequent classification task.

### Classification

The derived drug-consumption classes are related to personality characteristics through multinomial classification.

A cross-validated multinomial LASSO model is used to perform regularization and variable selection while reducing the risk of overfitting. Predictive performance is evaluated using a held-out test set and a confusion matrix.

The final model achieves approximately **74.0% classification accuracy**, indicating that personality characteristics contain useful predictive information about the derived drug-consumption profiles.

The analysis therefore forms an integrated statistical-learning pipeline:

**personality profiling --> dimensionality reduction of drug consumption --> weighted clustering --> response construction --> regularized multiclass prediction**

---

## Summary

**Statistical Learning – Drug Consumption & Personality:** Identified personality profiles using hierarchical and k-means clustering, constructed interpretable drug-consumption classes through PCA and Varimax rotation (`psych`, `nFactors`, `EFA.dimensions`), dependency-weighted Manhattan distance and PAM clustering (`cluster`, `factoextra`), and predicted these classes using cross-validated multinomial LASSO (`glmnet`), achieving approximately 74% classification accuracy.
