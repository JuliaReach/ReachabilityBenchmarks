
## declare an array variable
declare -a arr=("01" "005" "001" "0005")
for i in "${arr[@]}"
do
  echo "*************************$i***************************************"
  spaceex -g filtered_oscillator.64.$i.cfg -m filtered_oscillator.xml -v l -o out.gen
done