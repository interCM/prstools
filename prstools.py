import os
import numpy as np
import pandas as pd

DOSAGES = np.array([2, -1, 1, 0])


def file_len(fname):
    """
    Returns a number of lines in the file.
    """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def check_valid_bed(bed_fname):
    """
    Checks whether bed_fname is a valid bed file.
    See https://www.cog-genomics.org/plink2/formats#bed for bed format specs.
    """
    b = open(bed_fname, 'rb').read(3)
    if b != bytes([0x6c, 0x1b, 0x01]):
        raise RuntimeError(f"{bed_fname} file is not a valid bed file!")


def read_bed(bed_fname, n_snps, n_samples):
    """
    bed_fname - a name of a bed file.
    """
    n_blocks = n_snps
    block_len = n_samples//4 + (1 if n_samples%4 else 0)
    buf = open(bed_fname, 'rb')
    buf.read(3) # skip 'magic' bytes

    for _ in range(n_blocks):
        a1_dosages = _read_block(buf, block_len)
        a1_dosages = a1_dosages[:n_samples] # drop auxiliary format tail
        yield a1_dosages



def _read_block(buf, block_len):
    """
    Returns dosages (array) of 1st allele in bim file for the current sample.
    Contains elements {2, -1, 1, 0},
    -1 - for missing genotypes.
    """
    block = buf.read(block_len)
    ind = [(b >> 2*i) & 3 for b in block for i in range(4)]
    return DOSAGES[ind]



def estimate_prs(plink_root_fname, sumstats_fname, p_threshold=1):
    """
    """
    print(f"Reading {sumstats_fname} sumstats file")
    ss_df = pd.read_table(sumstats_fname, delim_whitespace=True)
    print(f"{len(ss_df)} SNPs in total")
    ss_df = ss_df[ss_df.P <= p_threshold]
    print(f"{len(ss_df)} SNPs with p-values <= {p_threshold}")

    print(f"Reading {plink_root_fname} plink files")
    bed_fname = f"{plink_root_fname}.bed"
    fam_fname = f"{plink_root_fname}.fam"
    bim_fname = f"{plink_root_fname}.bim"
    check_valid_bed(bed_fname)


    snps = pd.read_table(bim_fname, squeeze=True, header=None, delim_whitespace=True,
        names=["CHR", "SNP", "POS", "BP", "A1", "A2"], usecols=["SNP"])
    n_snps = len(snps)
    print(f"{n_snps} SNPs in bim file")

    fam_df = pd.read_table(fam_fname, header=None, delim_whitespace=True, 
        names=["FID", "IID", "FATHER", "MOTHER", "SEX", "PHEN"],
        usecols=["FID", "PHEN"])
    n_samples = len(fam_df)
    print(f"{n_samples} sampels in fam file")
    fam_df["PRS"] = 0

    ss_df = ss_df[ss_df.SNP.isin(snps)]
    print(f"{len(ss_df)} SNPs from sumstats file are in plink bim file")
    used_snps = frozenset(ss_df.SNP)

    effect_dict = dict(zip(ss_df.SNP, ss_df.BETA))
    i = 0
    for snp, dosages in zip(snps, read_bed(bed_fname, n_snps, n_samples)):
        if snp in used_snps:
            dosages[dosages == -1] = 0 # ignore effect of missing genotypes
            effect = effect_dict[snp]
            fam_df["PRS"] += effect*dosages
            i += 1
        # if i ==10:
            # break
    print(i)
    print(fam_df.head())





if __name__ == "__main__":
    plink_root_fname = "PRSice_toy_data/TOY_TARGET_DATA"
    sumstats_fname = "PRSice_toy_data/TOY_BASE_GWAS_BETA.assoc"
    estimate_prs(plink_root_fname, sumstats_fname)
