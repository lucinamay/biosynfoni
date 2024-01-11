from collections import Counter

import numpy as np
from scipy.stats import nbinom
import matplotlib.pyplot as plt


def check_cooccurrence(array):
    # Count the occurrences of each cell
    counts = Counter(array)

    # Calculate the probabilities of each cell
    total = len(array)
    probabilities = {cell: count / total for cell, count in counts.items()}

    # Calculate the joint probabilities
    joint_probabilities = Counter((array[i], array[i + 1]) for i in range(total - 1))

    # Check for co-occurrence
    co_occurrences = {}
    for (cell1, cell2), joint_count in joint_probabilities.items():
        joint_probability = joint_count / (total - 1)
        individual_product = probabilities[cell1] * probabilities[cell2]
        if joint_probability > individual_product:
            co_occurrences[(cell1, cell2)] = joint_probability / individual_product

    return co_occurrences


# # ============== Kolmogorov-Smirnov Test ==============================
# # test if the counts follow a negative binomial distribution
# from scipy.stats import kstest

# # Generate some data
# data = np.random.negative_binomial(1, 0.1, 1000)

# # Fit the negative binomial distribution
# params = nbinom.fit(data)

# # Perform the Kolmogorov-Smirnov test
# D, p_value = kstest(data, 'nbinom', args=params)

# print(f'D = {D}, p-value = {p_value}')


# ============== visualize the co-occurrence ==============================
# import seaborn as sns
# import pandas as pd
# import matplotlib.pyplot as plt

# # Assuming co_occurrences is your dictionary of co-occurring cells
# co_occurrences = check_cooccurrence(array)

# # Convert the dictionary to a DataFrame
# df = pd.DataFrame.from_dict(co_occurrences, orient='index', columns=['Co-occurrence'])

# # Create a pivot table for the heatmap
# pivot = df.reset_index().pivot('level_0', 'level_1', 'Co-occurrence')

# # Draw the heatmap
# sns.heatmap(pivot, cmap='YlGnBu')
# plt.title('Co-occurrence of Cells')
# plt.show()

# ======= Chi-Square Test =======
# from scipy.stats import chi2_contingency

# # Create a contingency table
# contingency_table = pd.crosstab(array[:-1], array[1:])

# # Perform the Chi-Square test
# chi2, p_value, _, _ = chi2_contingency(contingency_table)

# print(f'Chi-square = {chi2}, p-value = {p_value}')

# ===========================================================


def main():
    # data = np.loadtxt("data.txt", delimiter=",", skiprows=1)

    # Test the function
    # array = np.random.negative_binomial(1, 0.1, 1000)

    print(check_cooccurrence(array))
    # plot the array
    plt.plot(array)
    plt.show()
    # print the co-occurrence


if __name__ == "__main__":
    main()
