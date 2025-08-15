
import pandas as pd
import sys

inputfile=sys.argv[1]
outputfile=sys.argv[2]
# Load your summary table
df = pd.read_csv(inputfile)  # Replace with your actual file

# Define categories (columns)
categories = ["summary_lrhq", "summary_lrhqae", "summary_mm_ont", "summary_mm_sr", "summary_wm"]

# Define the condition: "best" or "best;other" in all categories
intersection_condition = df[categories].applymap(lambda x: x in ["best", "best;other"]).all(axis=1)

# Extract the matching read IDs
matching_reads = df.loc[intersection_condition, "read_id"]

# Save to a .txt file
matching_reads.to_csv(sys.argv[2], index=False, header=False)

print(f"Saved {len(matching_reads)} read IDs to best_best_other_reads.txt")
