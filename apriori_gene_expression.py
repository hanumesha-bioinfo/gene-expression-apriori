from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori, association_rules
import pandas as pd


# Gene expression transactions
transactions = [
    ['A', 'B', 'C'],
    ['A', 'B'],
    ['A', 'C'],
    ['B', 'C'],
    ['A', 'B', 'C'],
    ['A', 'B']
]


# Convert transactions into one-hot encoded format
te = TransactionEncoder()
te_data = te.fit(transactions).transform(transactions)

df = pd.DataFrame(te_data, columns=te.columns_)

print("Transaction Data:")
print(df)


# Apply Apriori algorithm with minimum support of 50%
frequent_itemsets = apriori(
    df,
    min_support=0.50,
    use_colnames=True
)

print("\nFrequent Itemsets:")
print(frequent_itemsets)


# Generate association rules with minimum confidence of 70%
rules = association_rules(
    frequent_itemsets,
    metric="confidence",
    min_threshold=0.70
)


print("\nAssociation Rules:")

if not rules.empty:

    print(
        rules[
            [
                "antecedents",
                "consequents",
                "support",
                "confidence",
                "lift"
            ]
        ]
    )


    # Identify strongest gene association
    strongest_rule = rules.loc[rules["lift"].idxmax()]

    print("\nStrongest Gene Association:")
    print(strongest_rule)

else:
    print("No association rules found.")
