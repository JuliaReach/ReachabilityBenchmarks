cd ./small/
for (( model_iter=2; model_iter<128; model_iter*=2 ))
do
  echo "*************************$model_iter***************************************"
  sspaceex -g filtered_oscillator.$model_iter.cfg -m filtered_oscillator.xml -v l -o out.gen
done
cd ../
cd ./128/
echo "*************************128***************************************"
sspaceex -g filtered_oscillator.128.cfg -m filtered_oscillator_128.xml -v l -o out.gen
cd ../
cd ./196/
echo "*************************196***************************************"
sspaceex -g filtered_oscillator.196.cfg -m filtered_oscillator_196.xml -v l -o out.gen
cd ../
cd ./256/
echo "*************************256***************************************"
sspaceex -g filtered_oscillator.256.cfg -m filtered_oscillator_256.xml -v l -o out.gen
cd ../
