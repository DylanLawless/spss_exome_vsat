#!/usr/bin/env bash
set -euo pipefail

bcftools="$HOME/web/tools/bcftools/bcftools"
post_dir="${1:-$HOME/web/spss_exome_vsat/data/post_ppi}"
out="${2:-$HOME/web/spss_exome_vsat/data/variant_qual}"

mkdir -p "$out"
[[ -x "$bcftools" ]] || { echo "bcftools not executable: $bcftools" >&2; exit 1; }

process_vcf() {
    local source="$1" vcf="$2" symbol_index kir_count

    [[ -f "$vcf" ]] || { echo "VCF not found: $vcf" >&2; exit 1; }

    # Derive the SYMBOL column from this VCF's VEP CSQ header.
    symbol_index=$(
        "$bcftools" view -h "$vcf" |
            awk '/^##INFO=<ID=CSQ,/ {
                line = $0
                sub(/^.*Format: /, "", line)
                n = split(line, f, "|")
                for (i = 1; i <= n; i++) {
                    gsub(/[^[:alnum:]_].*$/, "", f[i])
                    if (f[i] == "SYMBOL") { print i; exit }
                }
            }'
    )

    [[ "$symbol_index" =~ ^[0-9]+$ ]] || {
        echo "Missing CSQ SYMBOL annotation in $vcf" >&2
        exit 1
    }

    # Site-level quality, retaining annotated KIR genes.
    {
        printf 'CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tMQ\tQD\tFS\tSOR\tVQSLOD\tDP\tAC\tAF\tAN\tKIR_GENES\tIS_KIR\n'

        "$bcftools" view -m2 -M2 -v snps -Ou "$vcf" |
            "$bcftools" query -u -f \
'%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/MQ\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/VQSLOD\t%INFO/DP\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/CSQ\n' |
            awk -F '\t' -v OFS='\t' -v idx="$symbol_index" '
            {
                genes = ""
                n = split($16, annotations, ",")
                for (i = 1; i <= n; i++) {
                    split(annotations[i], fields, "|")
                    gene = fields[idx]

                    # Named KIR family genes, excluding KIR pseudogenes.
                    if (gene ~ /^KIR[0-9]/ &&
                        gene !~ /^KIR[0-9]DP/ &&
                        seen[gene] != NR) {
                        seen[gene] = NR
                        genes = genes (genes == "" ? "" : ";") gene
                    }
                }
                $16 = (genes == "" ? "." : genes)
                print $0, (genes == "" ? 0 : 1)
            }'
    } > "$out/site_quality_${source}.tsv"

    # Carrier-level quality; no rows for reference-only genotypes.
    {
        printf 'CHROM\tPOS\tREF\tALT\tSAMPLE\tGT\tDP\tGQ\tAD_REF\tAD_ALT\tAB\n'

        "$bcftools" view -m2 -M2 -v snps -Ou "$vcf" |
            "$bcftools" query -u -f \
'[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%DP\t%GQ\t%AD\n]' |
            awk -F '\t' -v OFS='\t' '
            {
                n = split($6, alleles, /[\/|]/)
                is_carrier = 0
                for (j = 1; j <= n; j++) {
                    if (alleles[j] == "1") is_carrier = 1
                }
                if (!is_carrier) next

                split($9, ad, ",")
                ref = (ad[1] ~ /^[0-9]+$/) ? ad[1] : "."
                alt = (ad[2] ~ /^[0-9]+$/) ? ad[2] : "."
                ab = "."
                if (ref != "." && alt != "." && ref + alt > 0)
                    ab = alt / (ref + alt)

                print $1, $2, $3, $4, $5, $6, $7, $8, ref, alt, ab
            }'
    } > "$out/carrier_quality_${source}.tsv"

    kir_count=$(awk -F '\t' 'NR > 1 && $17 == 1 { n++ } END { print n + 0 }' \
        "$out/site_quality_${source}.tsv")
    printf '%-12s %s KIR SNV records\n' "$source" "$kir_count"
}

for module in 22 586 836; do
    process_vcf "MCL_${module}" \
        "$post_dir/bcftools_gatk_norm_maf01.recode_vep_conda_impact_MCL_${module}.vcf.gz"
done

printf 'Outputs: %s\n' "$out"

