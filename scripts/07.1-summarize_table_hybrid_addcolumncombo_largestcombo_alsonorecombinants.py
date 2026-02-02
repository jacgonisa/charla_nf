from collections import Counter
import sys

def parse_crossovers(file_path, output_file):
    """
    Parse crossover data and add summary column.
    Optimized for speed with buffered writing and minimal string operations.
    """
    crossover_counts = Counter()

    # Use buffered writing for better performance
    with open(file_path, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            columns = line.rstrip('\n').split('\t')

            # Make sure there are at least 3 columns (we need columns[2])
            if len(columns) < 3:
                f_out.write(line.rstrip('\n') + '\tNA\n')
                continue

            # Get the crossover information from the third column
            crossovers = columns[2].split(',')

            # Extract the unique regions (the last 3 characters of each entry)
            # Using set comprehension for speed
            regions = sorted({x[-3:] for x in crossovers if x})

            if len(regions) > 1:
                # If there is recombination (more than one region), join them with a dash.
                abbreviation_combination = "-".join(regions)
                crossover_counts[abbreviation_combination] += 1
            elif len(regions) == 1:
                # If there's only one region (i.e. no recombination), use that region's abbreviation.
                abbreviation_combination = regions[0]
                crossover_counts[abbreviation_combination] += 1
            else:
                abbreviation_combination = "NA"

            # Write directly to file instead of building list in memory
            f_out.write(line.rstrip('\n') + f'\t{abbreviation_combination}\n')

    return crossover_counts

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    output_file = sys.argv[2]
    crossover_frequencies = parse_crossovers(file_path, output_file)

    # Print the combination counts
    for combo, freq in crossover_frequencies.most_common():
        print(f"{combo}: {freq}")

