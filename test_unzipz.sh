for i in ${1}/*.gz; do
  old_name=$(basename ${i} | rev | cut -d'.' -f3- | rev)
  new_name=$(echo ${old_name} | tr -d '[],').fasta
  dir_name=$(dirname ${i})
  echo -e "${i}\n${old_name}\n${new_name}\n${dir_name}"
  gunzip ${i}
  tax_genus=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f2 | tr -d '[],')
  tax_species=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f3 | tr -d '[],')
  echo "Taxes: ${tax_genus}:${tax_species}"
  mv ${dir_name}/${old_name} ${dir_name}/${tax_genus}_${tax_species}_${new_name}
done < "${1}"
