for i in "${1}/ANI/localANIDB_REFSEQ/*.gz"; do
  old_name=$(basename ${i} | cut -d'.' -f1,2)
  new_name=$(echo ${old_name} | tr -d '[],')
  dir_name=$(dirname ${i})
  echo "${old_name}-${new_name}-${dirname}"
  echo "unzipping ${i}"
  gunzip ${i}
  tax_genus=$(head -n1 "${dir_name}/${old_name}.fna" | cut -d' ' -f2 | tr -d '[],')
  tax_species=$(head -n1 "${dir_name}/${old_name}.fna" | cut -d' ' -f3 | tr -d '[],')
  echo "Taxes: ${tax_genus}:${tax_species}"
  echo "Moving ${dir_name}/${old_name}.fna to ${dir_name}/${tax_genus}_${tax_species}_${new_name}.fasta"
  mv ${dir_name}/${old_name}.fna ${dir_name}/${tax_genus}_${tax_species}_${new_name}.fasta
done
