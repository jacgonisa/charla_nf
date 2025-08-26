## ${projectDir}/scripts/08-plot_hybrid_summary.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# --- Argument Handling ---
# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python 08-plot_hybrid_summary.py <input_table.tsv> <output_prefix>")
    sys.exit(1)

input_file = sys.argv[1]
output_prefix = sys.argv[2]

# --- Data Loading and Processing ---
df = pd.read_csv(input_file, sep="\t")

# 1. Plot recombination events across thresholds
recomb_counts = df.groupby('threshold')['hybrid'].apply(lambda x: (x == "yes|2").sum()).reset_index()
recomb_counts.columns = ['threshold', 'recombination_events']

plt.figure(figsize=(8, 5))
sns.barplot(data=recomb_counts, x="threshold", y="recombination_events", palette="viridis")
plt.xlabel("Threshold")
plt.ylabel("Number of Recombination Events")
plt.title("Recombination Events Across Different Thresholds")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{output_prefix}_recombination_events.png", dpi=300)
plt.close() # Close the plot to free memory

# Save summary statistics to a file
with open(f"{output_prefix}_summary_stats.txt", "w") as f:
    f.write("Recombination Events Summary:\n")
    f.write(recomb_counts.describe().to_string())

# --- Filter data for subsequent plots ---
df_yes2 = df[df["hybrid"] == "yes|2"].copy() # Use .copy() to avoid SettingWithCopyWarning

# 2. Plot stacked bar chart of hybrid types
hybrid_counts_yes2 = df_yes2.groupby(['threshold', 'hybrid_type']).size().unstack(fill_value=0)

hybrid_counts_yes2.plot(kind='bar', stacked=True, colormap="viridis", figsize=(10, 6))
plt.xlabel("Threshold")
plt.ylabel("Count of Hybrid Types")
plt.title("Stacked Bar Plot of Hybrid Types by Threshold (Recombinants only)")
plt.legend(title="Hybrid Type", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{output_prefix}_hybrid_types_stacked.png", dpi=300, bbox_inches='tight')
plt.close()

# 3. Plot stacked bar chart of ectopic categories
df_yes2["ectopic_simplified"] = df_yes2["ectopic"].apply(lambda x: "yes" if x.startswith("yes") else x)
ectopic_counts = df_yes2.groupby(['threshold', 'ectopic_simplified']).size().unstack(fill_value=0)

ectopic_counts.plot(kind='bar', stacked=True, colormap="magma", figsize=(10, 6))
plt.xlabel("Threshold")
plt.ylabel("Count of Ectopic Categories")
plt.title("Stacked Bar Plot of Ectopic Categories by Threshold (Recombinants only)")
plt.legend(title="Ectopic Category", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{output_prefix}_ectopic_categories_stacked.png", dpi=300, bbox_inches='tight')
plt.close()
