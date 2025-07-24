from collections import Counter
import sys

def parse_crossovers(file_path, output_file):
    crossover_counts = Counter()
    updated_lines = []
    
    with open(file_path, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            # Make sure there are at least 3 columns (we need columns[2])
            if len(columns) < 3:
                updated_lines.append(line.strip() + '\tNA')
                continue
            
            # Get the crossover information from the third column
            crossovers = columns[2].split(',')
            
            # Extract the unique regions (the last 3 characters of each entry)
            regions = sorted(set(x[-3:] for x in crossovers if x))
            
            if len(regions) > 1:
                # If there is recombination (more than one region), join them with a dash.
                abbreviation_combination = "-".join(regions)
                crossover_counts.update([abbreviation_combination])
            elif len(regions) == 1:
                # If there's only one region (i.e. no recombination), use that region's abbreviation.
                abbreviation_combination = regions[0]
                crossover_counts.update([abbreviation_combination])
            else:
                abbreviation_combination = "NA"
            
            updated_lines.append(line.strip() + f'\t{abbreviation_combination}')
    
    # Write the output file with the new column
    with open(output_file, 'w') as f:
        for updated_line in updated_lines:
            f.write(updated_line + '\n')
    
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

