from pathlib import Path
import sys

def convert_mat_to_homer(input_file, fhandle):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Initialize variables
    motif_id = ""
    consensus = ""
    pwm = []
    header_found = False

    # Parse the .mat file
    for line in lines:
        line = line.strip()
        if line.startswith("ID"):
            motif_id = line.split("\t")[1]
        elif line.startswith("DE"):
            consensus = line.split("Consensus=")[1].split(",")[0]
        elif line.startswith("P0"):
            header_found = True
        elif header_found and line.startswith("XX"):
            break  # End of PWM section
        elif header_found:
            pwm.append(line.split("\t")[1:-1])  # Extract PWM values

    # Write header
    fhandle.write(f">{consensus}\t{motif_id}\t0\n") # make score 0, get continuous values
    # Write PWM
    for row in pwm:
        fhandle.write("\t".join(row) + "\n")

if __name__ == '__main__':
    mcross_mat_folder = Path(sys.argv[1])
    output_file = sys.argv[2]
    # Write to HOMER motif file
    with open(output_file, 'w') as fhandle:
        for input_file in mcross_mat_folder.glob('*.mat'):
            convert_mat_to_homer(input_file, fhandle)
        
