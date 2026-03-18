import argparse
import re
from collections import defaultdict

def main():
    # Bash betiğinden gelen argümanları yakala
    parser = argparse.ArgumentParser(description="Extract introns from GTF")
    parser.add_argument('--gtf', required=True, help="Path to input GTF file")
    parser.add_argument('--out', required=True, help="Path to output BED file")
    args = parser.parse_args()

    transcripts = {}
    exons = defaultdict(list)

    print(f"GTF okunuyor: {args.gtf}")
    with open(args.gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature = fields[2]
            if feature not in ('transcript', 'exon'):
                continue
            
            attrs = fields[8]
            
            # Sadece protein_coding genleri al
            if 'gene_type "protein_coding"' not in attrs:
                continue
            
            tid_m = re.search(r'transcript_id "([^"]+)"', attrs)
            gn_m  = re.search(r'gene_name "([^"]+)"', attrs)
            if not tid_m or not gn_m:
                continue
            
            tid   = tid_m.group(1)
            gene  = gn_m.group(1)
            chrom = fields[0]
            start = int(fields[3])
            end   = int(fields[4])
            strand = fields[6]
            
            if feature == 'transcript':
                transcripts[tid] = (chrom, gene, strand)
            elif feature == 'exon':
                exons[tid].append((start, end))

    print(f"Transcript sayısı: {len(transcripts)}")

    intron_lines = []
    for tid, (chrom, gene, strand) in transcripts.items():
        ex_list = sorted(exons.get(tid, []))
        
        if len(ex_list) < 2:
            continue
        
        for i in range(len(ex_list) - 1):
            # BED 0-based half-open format düzeltmesi
            bed_start = ex_list[i][1]
            bed_end = ex_list[i+1][0] - 1
            
            if bed_end > bed_start:
                # 6 sütunlu standart BED yapısı eklendi (skor olarak 0)
                intron_lines.append(f"{chrom}\t{bed_start}\t{bed_end}\t{gene}\t0\t{strand}\n")

    print("İntronlar sıralanıyor...")
    # Çıktıyı BEDTools'un beklediği gibi sırala
    intron_lines.sort(key=lambda x: (x.split('\t')[0], int(x.split('\t')[1])))

    # Kaydet
    with open(args.out, 'w') as f:
        f.writelines(intron_lines)

    print(f"Başarıyla kaydedildi: {args.out}")

if __name__ == "__main__":
    main()