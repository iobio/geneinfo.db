mane_gff=$1
grep tag=MANE_Select $mane_gff | grep '\ttranscript\t' | sed -rn 's/.*transcript_id=([A-Z.0-9]+).*gene_name=([A-Za-z.0-9]+).*Dbxref=RefSeq:([A-Z_0-9.]+).*/\2,\1,\3/p'
