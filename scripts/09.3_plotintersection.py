#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Plot intersections of read classifications across mapping methods using a parallel categories diagram."
    )
    parser.add_argument("csvfile", help="Input CSV file (e.g. summary_comparison.csv)")
    parser.add_argument("-o", "--outfile", default="parallel_categories.html", help="Output HTML file for the plot")
    args = parser.parse_args()

    # Read the CSV file (read_id is the index)
    df = pd.read_csv(args.csvfile, index_col="read_id")

    # Check if DataFrame is empty or has no columns
    if df.empty or len(df.columns) == 0:
        print("Warning: Input CSV is empty or has no data columns. Creating placeholder HTML.")
        # Create a simple HTML file with a message
        with open(args.outfile, 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head><title>No Data</title></head>
<body>
<h1>No Classification Data Available</h1>
<p>The summary table is empty. No reads were classified in this analysis.</p>
</body>
</html>""")
        print(f"Placeholder plot saved to {args.outfile}")
        return

    # Optionally, if cells contain multiple classifications (separated by semicolons),
    # we take only the first classification.
    for col in df.columns:
        df[col] = df[col].apply(lambda x: x.split(";")[0] if isinstance(x, str) and ";" in x else x)

    # Create a parallel categories plot.
    fig = px.parallel_categories(
        df,
        dimensions=df.columns.tolist(),
        color=df[df.columns[0]].astype("category").cat.codes,  # color by the first method's classification (as an example)
        color_continuous_scale=px.colors.sequential.Viridis,
        labels={col: col for col in df.columns},
        title="Intersection of Read Classifications Across Mapping Methods"
    )

    fig.write_html(args.outfile)
    print(f"Plot saved to {args.outfile}")

if __name__ == "__main__":
    main()

