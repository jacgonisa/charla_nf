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

