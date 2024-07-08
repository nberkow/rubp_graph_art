import sys
import os

export_dir = "../hit_list_exports/"

organism_set = set()
for res_file in os.listdir(export_dir):
    with open(f"{export_dir}{res_file}") as rf:
        for gene_line in rf:
            fields = gene_line.rstrip().split("\t")
            if len(fields) == 3 and "Organism" not in fields[0]:
                organism = fields[0]
                if not organism in organism_set:
                    transcript_id = fields[1].split("=")[-1]
                    print(f"{organism}\t{transcript_id}")
                    organism_set.add(organism)