<div align="center">

# Pairwise Independence Analysis of Data Quality Dimensions

*Supplementary data and analysis code for a manuscript under review at* **[ACM JDIQ](https://dl.acm.org/journal/jdiq)**

[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-3776ab?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
![anonymized for peer review](https://img.shields.io/badge/anonymized-for_peer_review-6b7280?style=flat-square)
[![license: CC BY 4.0](https://img.shields.io/badge/license-CC_BY_4.0-0085ca?style=flat-square)](LICENSE)

</div>

This repository contains the data artifacts, the raw contingency tables and the statistical analysis underlying the pairwise independence analysis reported in Section 5.3 and Appendix C of the manuscript.

## Data

The analysis rests on a catalog of 51 real-world data quality problems from the cultural heritage domain.
Each problem is annotated with at least one quality dimension, resulting in 185 assignments across 20/23 dimensions and thus C(20,2) = 190 pairwise comparisons.

| File | Description |
| --- | --- |
| [`problem-dimension-assignments.csv`](problem-dimension-assignments.csv) | Input data: one row per problem, with the primary and the further affected dimensions. Corresponds to Table 8 of the manuscript. |
| [`pairwise-independence-results.csv`](pairwise-independence-results.csv) | Results: one row per dimension pair, sorted by descending absolute phi. |
| [`sensitivity-analysis.txt`](sensitivity-analysis.txt) | Console output of the analysis, including the comparison between the two test strategies. |
| [`analyze_independence.py`](analyze_independence.py) | Analysis script producing the results from the input data. |

## Result columns

Each row of `pairwise-independence-results.csv` describes one pair of dimensions.

| Column | Description |
| --- | --- |
| `dimension_a`, `dimension_b` | The two dimensions of the pair |
| `freq_a`, `freq_b` | Marginal frequencies: problems affecting each dimension |
| `n11`, `n10`, `n01`, `n00` | Observed cell counts of the 2x2 contingency table. `n11` counts the problems affecting both dimensions, `n00` those affecting neither. |
| `exp_11`, `exp_10`, `exp_01`, `exp_00` | Expected cell counts under the null hypothesis of independence |
| `min_expected` | Smallest of the four expected cell counts. The chi-square approximation is considered admissible at `min_expected` >= 5, which holds for 14 of the 190 pairs. |
| `chi2`, `chi2_p` | Pearson chi-square statistic and its p-value, retained for comparison only |
| `fisher_p` | Two-sided p-value of Fisher's exact test, applied uniformly to all 190 pairs and the basis of the analysis |
| `fisher_p_holm` | Holm-adjusted p-value over all 190 comparisons, controlling the familywise error rate |
| `fisher_p_bh` | Benjamini-Hochberg adjusted p-value (q-value), controlling the false discovery rate |
| `phi` | Phi coefficient as a measure of effect size, ranging from -1 to +1 |

## Running the analysis

Using [uv](https://docs.astral.sh/uv/):

```bash
uv run analyze_independence.py
```

Using Python and pip:

```bash
pip install -r requirements.txt
python3 analyze_independence.py
```

## License

The contents of this repository are licensed under
[CC BY 4.0](LICENSE).
