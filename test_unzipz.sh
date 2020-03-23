for i in ${1}/*.gz; do
  old_name=$(basename ${i} | cut -d'.' -f1)
  new_name=$(echo ${old_name} | tr -d '[],')
  dir_name=$(dirname ${i})
  echo -e "${i}\n${old_name}.fna\n${new_name}.fasta\n${dir_name}"
  gunzip ${i}
  tax_genus=$(head -n1 "${dir_name}/${old_name}.fna" | cut -d' ' -f2 | tr -d '[],')
  tax_species=$(head -n1 "${dir_name}/${old_name}.fna" | cut -d' ' -f3 | tr -d '[],')
  echo "Taxes: ${tax_genus}:${tax_species}"
  mv ${dir_name}/${old_name}.fna ${dir_name}/${tax_genus}_${tax_species}_${new_name}.fasta
done < "${1}"
