for i in /raid5/MiSeqAnalysisFiles/test51_kp/1804718/ANI/localANIDB_REFSEQ/*.gz; do
  old_name=$(basename ${i} | rev | cut -d'.' -f2- | rev)
  new_name=$(echo ${old_name} | tr -d '[],')
  dir_name=$(dirname ${i})
  echo -e "${i}\n${old_name}\n${new_name}\n${dir_name}"
  gunzip ${i}
  tax_genus=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f2 | tr -d '[],')
  tax_species=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f3 | tr -d '[],')
  echo "Taxes: ${tax_genus}:${tax_species}"
  mv ${dir_name}/${old_name}.fna ${dir_name}/${tax_genus}_${tax_species}_${new_name}.fasta
done < "${1}"
