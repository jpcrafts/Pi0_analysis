import csv
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Convert Offsets_<run>_check.txt to a CSV file with block offsets and quality flags."
    )
    parser.add_argument("--run", type=int, required=True,
                        help="Run number to process (e.g. 123)")
    parser.add_argument("--input", type=str, default=None,
                        help="Optional input file name (default: Offsets_<run>_check.txt)")
    parser.add_argument("--output", type=str, default=None,
                        help="Optional output file name (default: Offsets_<run>.csv)")
    
    args = parser.parse_args()
    runnum = args.run
    input_filename = args.input if args.input else f"Offsets_{runnum}_check.txt"
    output_filename = args.output if args.output else f"Offsets_{runnum}.csv"
    
    offset_data = []  # List to store tuples of (offset_value, quality)

    try:
        with open(input_filename, 'r') as f:
            tokens = []
            for line in f:
                tokens.extend(line.strip().split())
    except FileNotFoundError:
        print(f"Error: File {input_filename} not found.")
        return

    # Process each token to extract the numeric value and quality flag.
    for token in tokens:
        quality = "OK"  # Default quality flag
        if token and token[-1] in [',', '!', '?']:
            if token[-1] == '!':
                quality = "FAIL"
            elif token[-1] == '?':
                quality = "LOW_STATS"
            token = token[:-1]  # Remove the trailing punctuation

        try:
            offset_value = float(token)
        except ValueError:
            continue  # Skip tokens that aren't valid numbers

        offset_data.append((offset_value, quality))

    if len(offset_data) != 1080:
        print(f"Warning: Expected 1080 offsets, but found {len(offset_data)}.")

    # Write the results to a CSV file.
    with open(output_filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Block", "Offset", "Quality"])
        for i, (offset, quality) in enumerate(offset_data):
            # Block numbering now starts at 0 (0 to 1079)
            csvwriter.writerow([i, offset, quality])

    print(f"CSV file '{output_filename}' has been created.")

if __name__ == '__main__':
    main()
