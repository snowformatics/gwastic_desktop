import subprocess


class GWAS:
    """
    GWAS class .
    """

    def vcf_to_bed(self, vcf_file, id_file, file_out, maf, geno):
        process = subprocess.Popen(["plink", "--vcf", vcf_file, "--make-bed", "--out", file_out,
                                    "--allow-extra-chr", "--set-missing-var-ids", "@:#", "--maf", maf,
                                    "--geno", geno, "--double-id"])#, "--keep", id_file])
        process.wait()