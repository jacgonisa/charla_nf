
import pandas as pd
import sys

inputfile=sys.argv[1]
outputfile=sys.argv[2]
# Load your summary table
df = pd.read_csv(inputfile)  # Replace with your actual file

# Define all possible categories (columns)
all_categories = ["summary_lrhq", "summary_lrhqae", "summary_mm_ont", "summary_mm_sr", "summary_wm"]

# Filter to only use categories that actually exist in the dataframe
categories = [cat for cat in all_categories if cat in df.columns]

if len(categories) == 0:
    print("Warning: No expected summary columns found in the input file.")
    print(f"Available columns: {list(df.columns)}")
    print("Creating empty output file.")
    # Create empty output file
    with open(outputfile, 'w') as f:
        pass
    sys.exit(0)

print(f"Using categories: {categories}")

# Define the condition: "best" or "best;other" in all available categories
intersection_condition = df[categories].map(lambda x: x in ["best", "best;other"]).all(axis=1)

# Extract the matching read IDs
matching_reads = df.loc[intersection_condition, "read_id"]

# Save to a .txt file
matching_reads.to_csv(sys.argv[2], index=False, header=False)

print(f"Saved {len(matching_reads)} read IDs to best_best_other_reads.txt")
