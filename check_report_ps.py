#!/usr/bin/env python3
import os
import re
import sys

def main(directory):
    # regex for lines like "Ps3_factor = -1"
    pat = re.compile(r'^Ps\d+_factor\s*=\s*(-?\d+)\s*$', re.IGNORECASE)

    # Allowed special PS factor values (1–7, plus 9 and 17)
    allowed_specials = set(range(1, 8)) | {9, 17, 129, 33}

    offending = []  # tuples: (filename, factors, unexpected, special_vals, special_count, neg_ones)
    processed = 0

    for fname in os.listdir(directory):
        path = os.path.join(directory, fname)
        if not os.path.isfile(path):
            continue

        # increment count and report every 1000 files
        processed += 1
        if processed % 1000 == 0:
            print(f"Processed {processed} files")

        factors = []
        with open(path, 'r') as f:
            for line in f:
                m = pat.match(line.strip())
                if m:
                    try:
                        factors.append(int(m.group(1)))
                    except ValueError:
                        factors.append(m.group(1))

        if not factors:
            # no PsN_factor lines → skip
            continue

        # Identify special values (1-7, 9, 17) and count -1's
        special_vals = [v for v in factors if v in allowed_specials]
        special_count = len(special_vals)
        neg_ones = factors.count(-1)

        # Any values not in allowed specials or -1
        unexpected = [v for v in factors if v not in allowed_specials and v != -1]

        # Check: exactly one special value, rest must be -1, and no unexpected
        if unexpected or special_count != 1 or neg_ones != len(factors) - 1:
            offending.append((fname, factors, unexpected, special_vals, special_count, neg_ones))

    # write detailed summary into cwd
    summary_path = os.path.join(os.getcwd(), 'daq_summary.txt')
    with open(summary_path, 'w') as out:
        out.write('DAQ Configuration errors\n')
        out.write('========================\n\n')
        for fname, factors, unexpected, special_vals, special_count, neg_ones in offending:
            out.write(f'File: {fname}\n')
            out.write(f'  All factors:        {factors}\n')
            if unexpected:
                out.write(f'  Unexpected vals:    {unexpected}\n')
            out.write(f'  Special values:     {special_vals}\n')
            out.write(f'  Special count:      {special_count}\n')
            out.write(f'  -1 count:           {neg_ones}\n\n')

    # print a simple list of offending files
    if offending:
        print("Offending files:")
        for fn, *_ in offending:
            print(f"  {fn}")
        print(f"\nSummary written to {summary_path}")
    else:
        print("All files clean (exactly one special PS factor and all others -1, with no unexpected values).")

if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    main(target_dir)
