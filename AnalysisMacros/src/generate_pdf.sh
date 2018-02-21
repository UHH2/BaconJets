file_name=controlplotsSiBCD_20Feb
if [ "$#" = "1" ]
   then file_name=controlplots${1}_trgM_31Jan
fi
rm ${file_name}.tex
py_out=$(./generate_tex_siL2Res.py "$file_name")
pdflatex ${file_name}.tex
rm ${file_name}.[tosnla][^e]*
echo -e "\nthe python output was:"
echo "$py_out"
