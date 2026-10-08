# Gene Expression Apriori Analysis

## Project Description

This project applies the Apriori algorithm to gene-expression transactions to identify frequent gene combinations and generate association rules.

The analysis uses a minimum support of 50% and a minimum confidence of 70%.

## Gene Expression Transactions

The dataset contains six gene-expression transactions:

- S1 = {A, B, C}
- S2 = {A, B}
- S3 = {A, C}
- S4 = {B, C}
- S5 = {A, B, C}
- S6 = {A, B}

## Methodology

The project performs the following steps:

1. Define gene-expression transactions.
2. Convert transactions into one-hot encoded format.
3. Apply the Apriori algorithm.
4. Identify frequent itemsets using 50% minimum support.
5. Generate association rules using 70% minimum confidence.
6. Calculate support, confidence, and lift.
7. Identify the strongest association based on the highest lift.

## Frequent Itemsets

### Frequent 1-itemsets

| Itemset | Support |
|---|---:|
| A | 83.33% |
| B | 83.33% |
| C | 66.67% |

### Frequent 2-itemsets

| Itemset | Support |
|---|---:|
| A, B | 66.67% |
| A, C | 50.00% |
| B, C | 50.00% |

### Frequent 3-itemsets

There are no frequent 3-itemsets because:

A, B, C occurs in 2 out of 6 transactions.

Support = 2/6 = 33.33%

This is below the minimum support of 50%.

## Association Rules

| Rule | Support | Confidence | Lift |
|---|---:|---:|---:|
| A → B | 66.67% | 80% | 0.96 |
| B → A | 66.67% | 80% | 0.96 |
| C → A | 50.00% | 75% | 0.90 |
| C → B | 50.00% | 75% | 0.90 |

## Strongest Gene Association

The strongest association based on the highest lift is between A and B.

A → B and B → A both have:

- Support = 66.67%
- Confidence = 80%
- Lift = 0.96

The two rules are tied for the highest lift.

Because the lift is slightly below 1, the data do not indicate a positive association between A and B compared with independence.

## Requirements

- Python 3
- pandas
- mlxtend

## Installation

Create and activate a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
# gene-expression-apriori
