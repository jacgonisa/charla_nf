import pandas as pd
import sys



myfile=sys.argv[1]
outfile=sys.argv[2]

# Read the data from a TSV file
# Adjust "your_file.tsv" to your actual filename.
# If your file already has headers, remove header=None and names.
df = pd.read_csv(myfile, sep="\t", header=None, names=["read", "start", "end", "label"])

# Sort the DataFrame by read and start coordinate
df = df.sort_values(by=["read", "start"])

results = []
# Group by each read
for read, group in df.groupby("read"):
    group = group.reset_index(drop=True)
    # For each consecutive pair, compute the width
    for i in range(len(group) - 1):
        gap = group.loc[i + 1, "start"] - group.loc[i, "end"]
        # If segments overlap, set gap to 0
        gap = gap if gap > 0 else 0
        segment_pair = f"{i + 1}_{i + 2}"
        results.append({"read": read, "segment_pair": segment_pair, "width": gap})

result_df = pd.DataFrame(results)


result_df.to_csv(outfile, sep ="\t", index=False)

