import pysam
import csv
import sys

ID=sys.argv[1]

# Input file paths
bam_path = "/n/data1/hms/dbmi/park/jiny/SMaHT/Tissue6/0.Data/All_LR_Tissue/"
bam_file = bam_path + ID  + "_all_pacbio.bam"  # Path to the BAM file
positions_file = sys.argv[2]  # Tab-delimited file with CHROM, Var_pos, Var_ref, Var_alt, Germ_pos, Germ_ref, Germ_alt

# Output file path
output_file = positions_file[:-4] + '_Readcnts.tsv'

# Open the BAM file
bam = pysam.AlignmentFile(bam_file, "rb")

# Read positions and analyze variants
with open(positions_file, "r") as pos_file, open(output_file, "w") as out_file:
    reader = csv.reader(pos_file, delimiter="\t")
    writer = csv.writer(out_file, delimiter="\t")

    # Write header to the output file
    header = next(reader) + ["total_read", "case_both", "case_var_only", "case_germ_only", "case_none"]
    writer.writerow(header)

    for line in reader:
        chrom, var_pos, var_ref, var_alt, germ_pos, germ_ref, germ_alt = (
            line[0],
            int(line[1]),
            line[2],
            line[3],
            line[4],
            line[5],
            line[6],
        )

        # Check if germ_pos is "NA"
        if germ_pos == "NA":
            # Write NA for all additional columns
            writer.writerow(line + ["NA", "NA", "NA", "NA", "NA"])
            continue

        # Ensure Var_pos and Germ_pos are integers
        germ_pos = int(germ_pos)
        start_pos = min(var_pos, germ_pos)
        end_pos = max(var_pos, germ_pos)

        # Initialize counters for each case
        case_counts = {"case_both": 0, "case_var_only": 0, "case_germ_only": 0, "case_none": 0}
        total_read = 0

        # Fetch reads covering both positions
        reads_var = {read.query_name: read for read in bam.fetch(chrom, var_pos - 1, var_pos)}
        reads_germ = {read.query_name: read for read in bam.fetch(chrom, germ_pos - 1, germ_pos)}

        # Intersect reads that span both positions
        common_reads = set(reads_var.keys()).intersection(reads_germ.keys())

        for read_name in common_reads:
            read = reads_var[read_name]  # Read from either reads_var or reads_germ (identical query_name)
            total_read += 1

            # Extract read sequence and reference positions
            read_seq = read.query_sequence
            read_pos = read.get_reference_positions(full_length=True)

            # Initialize flags for Var_pos and Germ_pos variants
            var_is_alt = False
            germ_is_alt = False

            # Check Var_pos
            if var_pos - 1 in read_pos:
                idx_var = read_pos.index(var_pos - 1)
                base_var = read_seq[idx_var]
                if base_var == var_alt:
                    var_is_alt = True

            # Check Germ_pos
            if germ_pos - 1 in read_pos:
                idx_germ = read_pos.index(germ_pos - 1)
                base_germ = read_seq[idx_germ]
                if base_germ == germ_alt:
                    germ_is_alt = True

            # Classify the read into a case
            if var_is_alt and germ_is_alt:
                case_counts["case_both"] += 1
            elif var_is_alt:
                case_counts["case_var_only"] += 1
            elif germ_is_alt:
                case_counts["case_germ_only"] += 1
            else:
                case_counts["case_none"] += 1

        # Write the result for the current line
        writer.writerow(
            line + [
                total_read,
                case_counts["case_both"],
                case_counts["case_var_only"],
                case_counts["case_germ_only"],
                case_counts["case_none"],
            ]
        )

# Close the BAM file
bam.close()

print(f"Results saved to {output_file}")
