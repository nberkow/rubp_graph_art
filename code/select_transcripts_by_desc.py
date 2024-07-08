from Bio import SeqIO


def get_best_desc(description_list):

    best_desc = { "rbcS" : ("",""),
                  "rbcL" : ("","")}
    
    best_score = {"rbcS" : 0,
                  "rbcL" : 0}
    
    best_str_len = {"rbcS" : 10**9,
                    "rbcL" : 10**9}

    base_keywords = "rubp", "rubisco", "ribulose", "bisphosphate", "carboxylase"
    s_keywords = "rbcS", "small"
    l_keywords = "rbcL", "large"

    for d in description_list:
        gid, desc = d

        kw = 0
        for k in base_keywords:
            if k in desc:
                kw += 1
        
        rbs_kw = kw
        for k in s_keywords:
            if k in desc:
                rbs_kw += 5
        for k in l_keywords:
            if k in desc:
                rbs_kw -= 5

        if rbs_kw >= best_score["rbcS"]:
            if len(desc) < best_str_len["rbcS"]:
                best_desc["rbcS"] = d
                best_str_len["rbcS"] = len(desc)

        rbl_kw = kw
        for k in l_keywords:
            if k in desc:
                rbl_kw += 5
        for k in s_keywords:
            if k in desc:
                rbl_kw -= 5

        if rbl_kw >= best_score["rbcL"]:
            if len(desc) < best_str_len["rbcL"]:
                best_desc["rbcL"] = d
                best_str_len["rbcL"] = len(desc)

    return(best_desc)


if __name__ == "__main__":

    transcript_fa = "../transcripts/transcripts.fa"
    seq_by_id = {}
    descriptions_by_organism = {}

    fasta_sequences = SeqIO.parse(transcript_fa,'fasta')
    for s in fasta_sequences:
        desc_elements = s.description.split("|")
        gene_id = desc_elements[0]
        organism = desc_elements[2]
        gene_desc = desc_elements[3]

        seq_by_id[gene_id] = s
        if not organism in descriptions_by_organism:
            descriptions_by_organism[organism] = []
        descriptions_by_organism[organism].append((gene_id, gene_desc))

    rbcS_sequences = []
    rbcL_sequences = []

    for og in descriptions_by_organism:

        best_desc = get_best_desc(descriptions_by_organism[og])

        if best_desc["rbcS"][0] != "":
            rbcS_sequences.append(seq_by_id[best_desc["rbcS"][0]])
        if best_desc["rbcL"][0] != "":
            rbcL_sequences.append(seq_by_id[best_desc["rbcL"][0]])

    with open("../transcripts/rbcS.fa", 'w') as rbcS_fasta:
        SeqIO.write(rbcS_sequences, rbcS_fasta, "fasta")

    with open("../transcripts/rbcL.fa", 'w') as rbcL_fasta:
        SeqIO.write(rbcL_sequences, rbcL_fasta, "fasta")
