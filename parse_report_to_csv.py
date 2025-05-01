#!/usr/bin/env python3

import os
import re
import sys
import csv

def parse_report_file(fname):
    """
    Parse a single .report file to extract:
      - bcm4a_charge_uc (BCM4A Charge in microCoulombs),
      - computer_live_time (percent, from Ps6/TRIG6),
      - tracking_eff (E SING FID TRACK EFFIC main value),
      - tracking_eff_err (E SING FID TRACK EFFIC uncertainty),
      - hod_eff (3_of_4 EFF).

    Returns a dict:
      {
        "bcm4a_charge_uc": float or None,
        "computer_live_time": float or None,
        "tracking_eff": float or None,
        "tracking_eff_err": float or None,
        "hod_eff": float or None
      }
    """

    results = {
        "bcm4a_charge_uc": None,
        "computer_live_time": None,
        "tracking_eff": None,
        "tracking_eff_err": None,
        "hod_eff": None
    }

    if not os.path.isfile(fname):
        print(f"[parse] File not found: {fname}")
        return results

    with open(fname, "r") as f:
        for line in f:
            # Example line: BCM4A Charge: 26042.240 uC
            if "BCM4A Charge:" in line:
                match = re.search(r"Charge:\s*([\d.]+)", line)
                if match:
                    results["bcm4a_charge_uc"] = float(match.group(1))

            # Example line: Pre-Scaled Ps6 HMS Computer Live Time : 99.9730 %
            # or HMS TRIG6 Computer Live Time : 99.97 %
            if ("Computer Live Time" in line) and ("Ps6" in line or "TRIG6" in line):
                match = re.search(r":\s*([\d.]+)\s*%", line)
                if match:
                    results["computer_live_time"] = float(match.group(1))

            # Example line: E SING FID TRACK EFFIC : 0.9954 +- 0.0002
            if "E SING FID TRACK EFFIC" in line:
                match = re.search(r":\s*([\d.]+)\s*\+-\s*([\d.]+)", line)
                if match:
                    results["tracking_eff"] = float(match.group(1))
                    results["tracking_eff_err"] = float(match.group(2))

            # Example line: 3_of_4 EFF : 0.999381
            if "3_of_4 EFF" in line:
                match = re.search(r":\s*([\d.]+)", line)
                if match:
                    results["hod_eff"] = float(match.group(1))

    return results


def parse_run_and_segment(fname):
    """
    Attempt to extract run number and segment from a filename like:
      coin_NPS_HMS_report_<run>_<seg>_1_-1.report

    Returns (run, seg) as integers, or (None, None) if not found.
    """
    base = os.path.basename(fname)
    # Example base name: coin_NPS_HMS_report_4084_0_1_-1.report
    # We'll do a regex capturing coin_NPS_HMS_report_(\d+)_(\d+)_1_-1\.report
    match = re.match(r"coin_NPS_HMS_report_(\d+)_(\d+)_1_-1\.report", base)
    if match:
        run_str = match.group(1)
        seg_str = match.group(2)
        return int(run_str), int(seg_str)
    else:
        # If the file name is in some other format, can't parse run/seg
        return None, None


def main():
    # If no arguments, print usage
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <report_file1> [report_file2 ...]")
        sys.exit(1)

    # We'll store each row of data in a list
    rows = []

    for fname in sys.argv[1:]:
        run, seg = parse_run_and_segment(fname)
        vals = parse_report_file(fname)
        row = {
            "run": run,
            "segment": seg,
            "bcm4a_charge_uc": vals["bcm4a_charge_uc"],
            "computer_live_time": vals["computer_live_time"],
            "tracking_eff": vals["tracking_eff"],
            "tracking_eff_err": vals["tracking_eff_err"],
            "hod_eff": vals["hod_eff"]
        }
        rows.append(row)

    # Write CSV to 'report_summary.csv' (change if you want a different filename)
    out_csv = "report_summary.csv"
    with open(out_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(["run", "segment", "bcm4a_charge_uc",
                         "computer_live_time", "tracking_eff",
                         "tracking_eff_err", "hod_eff"])
        # Write each row
        for r in rows:
            writer.writerow([
                r["run"],
                r["segment"],
                r["bcm4a_charge_uc"],
                r["computer_live_time"],
                r["tracking_eff"],
                r["tracking_eff_err"],
                r["hod_eff"]
            ])

    print(f"Wrote {len(rows)} rows to {out_csv}.")


if __name__ == "__main__":
    main()
